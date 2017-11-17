UMIgraphIndex <- function(UMIunique, UMIs, levenshtein_filt = 6, levenshtein_max = 12)
    #Filtered UMI's are compared by computing the levenshtein distance to each other UMI available in the list.
    #The levenshtein distance is used to measure the difference between two strings.
    #It is defined as the minimum number of single character changes (substitutions, insertions, deletions) required to transform one string into the other.
    #UMIlevenshtein returns a list containing an integer vector of levenshtein distances for each UMI.
    
    #UMIs_sort = UMIs_sort$UMIunique
    #UMIs = UMIs$UMI
{
    
    #Getting the levenshtein distance between every UMI.
    UMI_levenshtein.out <- .UMIlevenshtein(UMIunique, levenshtein_max = levenshtein_max)
    
    sortedUMI <- vector("list", length(UMI_levenshtein.out))
    for (i in 1:length(UMI_levenshtein.out)){
        if (length((which(UMI_levenshtein.out[[i]]<=levenshtein_filt&sort(table(UMIs), decreasing=TRUE)[[i]]>=2*sort(table(UMIs), decreasing=TRUE)-1)))!=0){
            sortedUMI[[i]] <- rbind(i, which(UMI_levenshtein.out[[i]]<=levenshtein_filt&sort(table(UMIs), decreasing=TRUE)[[i]]>=2*sort(table(UMIs), decreasing=TRUE)-1))
        }else{
            
        }
    }
    
    sortedUMI.bind <- do.call(cbind, sortedUMI)
    UMI_mat_prep <- matrix(sortedUMI.bind, nc = 2, byrow = TRUE)
    for(i in 1:nrow(UMI_mat_prep)){
        for(a in 1:nrow(UMI_mat_prep)){ 
            if (sum(UMI_mat_prep[i,]==UMI_mat_prep_rev[a,]) > 1){
                UMI_mat_prep[a,] <- c(0,0)
            }
        }
    }
    UMI_mat_prep <- UMI_mat_prep[-c(which(UMI_mat_prep[,1]==0)),]
    UMI_graph <- make_graph(as.vector(t(UMI_mat_prep)), n = length(sortedUMI), directed = FALSE)
    
    UMI_cluster <- .UMIcluster(UMIunique = UMIunique, UMI_graph = UMI_graph)
    
    UMI_index <- .UMIindex(UMI_cluster = UMI_cluster, UMIs = UMIs)

    
    return(list(UMIgraph = UMI_graph, UMIindex = UMI_index))
}


.UMIlevenshtein <- function(UMIunique, levenshtein_max = 6)
{
    UMI_levenshtein <- vector("list", length(UMIunique))
    
    for (i in 1:length(UMIunique)){
        UMI_levenshtein[[i]] <- stringdist(UMIunique[i], UMIunique)
        UMI_levenshtein[[i]][i] <- levenshtein_max
    } 
    
    return(UMI_levenshtein)
}


.UMIcluster <- function(UMIunique, UMI_graph)
    
    #UMI's are split corresponding to their graph node.
    
{
    UMI_clusters <- split(UMIs_sort, membership(components(UMI_graph)))
    
    return(UMI_clusters)
}

.UMIindex <- function(UMI_cluster, UMIs)
    
    #Initial UMI's (before sorting and filtering for duplicates) are indexed to receive an useful index for the corresponding reads (sequences).
    #UMIindex takes the initial corresponding UMI-read index into account and assigns each UMI node the index numbers of the initial UMI's that are part of the node by either being a duplicate or being connected to the id by their alignment.
    #This final list of index numbers makes it possible to subgroup the reads by their UMI's.
{
    UMI_index <- vector("list", length(UMI_cluster))
    for (i in 1:length(UMI_cluster)){
        len <- 1
        for (a in 1:length(UMI_cluster[[i]])){
            id <- which(UMIs==UMI_cluster[[i]][a])
            for (b in 1:length(id)){
                UMI_index[[i]][len] <- id[b] 
                len <- len+1
            }
        }
    }
    
    return(UMI_index)
}
