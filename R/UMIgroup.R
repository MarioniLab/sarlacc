


# Function for both UMIs.
# id.filt needs to be specified.
# umiGroup adds the levenshtein distances of both UMI's up and then proceeds as in the initial UMIgroup function.

umiGroup <- function(id.filt, UMI1, UMI2 = NULL, UMI_length = 12, Levensh_threshold = 2)
{
    if(!is.null(UMI2)){
        sort.dist1 <- .umiLevensh(id.filt, UMI1, UMI_length)
        sort.dist2 <- .umiLevensh(id.filt, UMI2, UMI_length)
        
        UMI.distance <- matrix(unlist(sort.dist2$UMI_levenshtein.out), length(sort.dist2$UMI_levenshtein.out))+matrix(unlist(sort.dist1$UMI_levenshtein.out), length(sort.dist1$UMI_levenshtein.out))
        Levensh_threshold <- Levensh_threshold*2
        
        sortedUMI <- vector("list", length = nrow(UMI.distance))
        for (i in 1:nrow(UMI.distance)){
            if (length((which(UMI.distance[i,]<=Levensh_threshold)))!=0){
                sortedUMI[[i]] <- rbind(i, which(UMI.distance[i,]<=Levensh_threshold))
            }else{
                sortedUMI
            }
        }
        
    }else{
        sort.dist1 <- .umiLevensh(id.filt, UMI1, UMI_length)
        UMI.distance <- sort.dist1$UMI_levenshtein.out
        
        sortedUMI <- vector("list", length(sort.dist1$UMI_levenshtein.out))
        for (i in 1:length(sort.dist1$UMI_levenshtein.out)){
            if (length((which(UMI.distance[[i]]<=Levensh_threshold)))!=0){
                sortedUMI[[i]] <- rbind(i, which(UMI.distance[[i]]<=Levensh_threshold))
            }else{
                sortedUMI
            }
        }
    }
    
    # Connected UMIs are shown in the matrix and used to make a graph.
    sortedUMI.bind <- do.call(cbind, sortedUMI)
    if(!is.null(sortedUMI.bind)){
        UMI_mat_prep <- matrix(sortedUMI.bind, nc = 2, byrow = TRUE)
        UMI_mat_prep_rev <- cbind(UMI_mat_prep[,2], UMI_mat_prep[,1])
        
        for(i in 1:nrow(UMI_mat_prep)){
            for(a in 1:nrow(UMI_mat_prep)){ 
                if (sum(UMI_mat_prep[i,]==UMI_mat_prep_rev[a,]) > 1){
                    UMI_mat_prep[a,] <- c(0,0)
                }
            }
        }
        UMI_mat_prep <- UMI_mat_prep[-c(which(UMI_mat_prep[,1]==0)),]
        UMI_graph <- make_graph(as.vector(t(UMI_mat_prep)), n = length(sortedUMI), directed = FALSE)
        
        # Indexing the input UMIs to be able to subset the corresponding reads.
        UMI_index <- membership(components(UMI_graph))
    }else{
        UMI_index <- NULL
        
    }
    
    return(UMI.id = UMI_index)
    
}


.umiLevensh <- function(id.filt, UMI, UMI_length)
{
    id.nfilt <- as.integer(id.filt)
    UMI_filt <- UMI[id.nfilt]
    UMI_filtDNA <- DNAStringSet(UMI_filt)
    UMI_levenshtein.out <- .UMIlevenshtein(UMI_filtDNA, levenshtein_max = UMI_length)
    
    return(list(UMI_levenshtein.out = UMI_levenshtein.out, UMI_filtDNA = UMI_filtDNA, UMI_filt = UMI_filt))
}

.addDistance <- function(UMI1_levenshtein.out, UMI2_levenshtein.out)
{
    UMI.distance <- vector("list", length(UMI1_levenshtein.out))
    for(i in 1:length(UMI1_levenshtein.out)){
        UMI.distance[[i]] <- UMI1_levenshtein.out[[i]] + UMI2_levenshtein.out[[i]]
    }
    
    return(UMI.distance)
}

# Calculates levenshtein distance between UMIs.
.UMIlevenshtein <- function(UMI_filtDNA, levenshtein_max = 6)
{
    UMI_levenshtein <- vector("list", length(UMI_filtDNA))
    
    for (i in 1:length(UMI_filtDNA)){
        UMI_levenshtein[[i]] <- stringdist(UMI_filtDNA[i], UMI_filtDNA)
        UMI_levenshtein[[i]][i] <- levenshtein_max
    } 
    
    return(UMI_levenshtein)
}
