


# Function for both UMIs.
# id.filt needs to be specified.
# umiGroup adds the levenshtein distances of both UMI's up and then proceeds as in the initial UMIgroup function.

umiGroup <- function(id.filt, UMI1, UMI2 = NULL, Levensh_threshold = 2)
{
    UMI_members = NULL
    for(i in 1:length(id.filt)){
        id.nfilt <- as.integer(id.filt[[i]])
        UMI1_filt <- UMI1[id.nfilt]
        UMI1_filtDNA <- DNAStringSet(UMI1_filt)
        
        if(!is.null(UMI2)){
            UMI2_filt <- UMI2[id.nfilt]
            UMI2_filtDNA <- DNAStringSet(UMI2_filt)
        }else{
            UMI2_filtDNA <- NULL
        }
        
        UMI.id <- .umiClust(UMI1 = UMI1_filtDNA, UMI2 = UMI2_filtDNA, Levensh_threshold = Levensh_threshold)
        if(!is.null(UMI.id)){
            UMI_members <- append(UMI_members, paste(i, UMI.id))
        }
    }
    
    return(UMI_members)
}

.umiClust <- function(UMI1, UMI2 = NULL, Levensh_threshold = 2)
{
    if(!is.null(UMI2)){
        
        UMI1_levenshtein.out <- .UMIlevenshtein(UMI_filtDNA = UMI1)
        UMI2_levenshtein.out <- .UMIlevenshtein(UMI_filtDNA = UMI2)
        
        UMI.distance <- matrix(unlist(UMI2_levenshtein.out), length(UMI2_levenshtein.out))+matrix(unlist(UMI1_levenshtein.out), length(UMI1_levenshtein.out))
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
        UMI1_levenshtein.out <- .UMIlevenshtein(UMI_filtDNA = UMI1)
        UMI.distance <- UMI1_levenshtein.out
        
        sortedUMI <- vector("list", length(UMI1_levenshtein.out))
        for (i in 1:length(UMI1_levenshtein.out)){
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

# Calculates levenshtein distance between UMIs.
.UMIlevenshtein <- function(UMI_filtDNA)
{
    UMI_levenshtein <- vector("list", length(UMI_filtDNA))
    
    for (i in 1:length(UMI_filtDNA)){
        UMI_levenshtein[[i]] <- stringdist(UMI_filtDNA[i], UMI_filtDNA)
        UMI_levenshtein[[i]][i] <- length(UMI_filtDNA)
    } 
    
    return(UMI_levenshtein)
}
