


# Function for both UMIs.
# id.filt needs to be specified.
# umiGroup adds the levenshtein distances of both UMI's up and then proceeds as in the initial UMIgroup function.

umiGroup1 <- function(id.filt, UMI1, UMI2 = NULL, UMI_length = 12, Levensh_threshold = 2)
{
    if(!is.null(UMI2)){
        sort.dist1 <- .umiLevensh(id.filt, UMI1, UMI_length)
        sort.dist2 <- .umiLevensh(id.filt, UMI2, UMI_length)
        
        UMI.distance <- .addDistance(sort.dist1$UMI_levenshtein.out, sort.dist2$UMI_levenshtein.out)
        Levensh_threshold <- Levensh_threshold*2
        sorted_UMI.filt <- sort(table(sort.dist1$UMI_filt), decreasing=TRUE) + sort(table(sort.dist2$UMI_filt), decreasing=TRUE)
        
    }else{
        sort.dist1 <- .umiLevensh(id.filt, UMI1, UMI_length)
        UMI.distance <- sort.dist1$UMI_levenshtein.out
        UMI_filt <- sort.dist1$UMI_filt
        sorted_UMI.filt <- sort(table(sort.dist1$UMI_filt), decreasing=TRUE)
    }
    
    sortedUMI <- vector("list", length(UMI.distance))
    substract <- sum(c(!is.null(UMI2), !is.null(UMI1)))
    for (i in 1:length(UMI.distance)){
        if (length((which(UMI.distance[[i]]<=Levensh_threshold&sorted_UMI.filt[[i]]>=2*sorted_UMI.filt-substract)))!=0){
            sortedUMI[[i]] <- rbind(i, which(UMI.distance[[i]]<=Levensh_threshold&sorted_UMI.filt[[i]]>=2*sorted_UMI.filt-substract))
        }else{
            sortedUMI
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
        if(!is.null(UMI2)){
            UMI1_index <- .umiIndex(UMI_sort = sort.dist1$UMI_sort, UMI_graph = UMI_graph, UMI_filt = sort.dist1$UMI_filt)
            UMI2_index <- .umiIndex(UMI_sort = sort.dist2$UMI_sort, UMI_graph = UMI_graph, UMI_filt = sort.dist2$UMI_filt)
            UMI_index <- list(UMI1_index, UMI2_index)
        }else{
            UMI_index <- .umiIndex(UMI_sort = sort.dist1$UMI_sort, UMI_graph = UMI_graph, UMI_filt = sort.dist1$UMI_filt)
        }
    }else{
        UMI_index <- NULL
        
    }
    
    return(UMI.id = UMI_index)
    
}



.umiIndex <- function(UMI_sort, UMI_graph, UMI_filt)
{
    UMI_cluster <- split(UMI_sort, membership(components(UMI_graph)))
    
    UMI_index <- vector("list", length(UMI_cluster))
    for (i in 1:length(UMI_cluster)){
        len <- 1
        for (a in 1:length(UMI_cluster[[i]])){
            id <- which(UMI_filt==UMI_cluster[[i]][a])
            for (b in 1:length(id)){
                UMI_index[[i]][len] <- id[b] 
                len <- len+1
            }
        }
    }
    return(UMI_index)
}

.umiLevensh <- function(id.filt, UMI, UMI_length)
{
    id.nfilt <- as.integer(id.filt)
    UMI_filt <- UMI[id.nfilt]
    UMI_sort <- DNAStringSet(names(sort(table(UMI_filt), decreasing=TRUE)))
    UMI_levenshtein.out <- .UMIlevenshtein(UMI_sort, levenshtein_max = UMI_length)
    
    return(list(UMI_levenshtein.out = UMI_levenshtein.out, UMI_sort = UMI_sort, UMI_filt = UMI_filt))
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
.UMIlevenshtein <- function(UMI_sort, levenshtein_max = 6)
{
    UMI_levenshtein <- vector("list", length(UMI_sort))
    
    for (i in 1:length(UMI_sort)){
        UMI_levenshtein[[i]] <- stringdist(UMI_sort[i], UMI_sort)
        UMI_levenshtein[[i]][i] <- levenshtein_max
    } 
    
    return(UMI_levenshtein)
}
