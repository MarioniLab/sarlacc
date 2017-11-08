UMIgraph <- function(UMI_unique_sort, levenshtein_filt = 1, UMIs)
    #Filtered UMI's are compared by computing the levenshtein distance to each other UMI available in the list.
    #The levenshtein distance is used to measure the difference between two strings.
    #It is defined as the minimum number of single character changes (substitutions, insertions, deletions) required to transform one string into the other.
    #UMIlevenshtein returns a list containing an integer vector of levenshtein distances for each UMI.
{
    
    #Getting the levenshtein distance between every UMI.
    UMI_levenshtein.out <- .UMIlevenshtein(UMI_unique_sort = UMI_unique_sort, levenshtein_max = 6)
    
    sortedUMI <- vector("list", length(UMI_levenshtein.out))
    for (i in 1:length(UMI_levenshtein.out)){
        if (length((which(UMI_levenshtein.out[[i]]<=levenshtein_filt&sort(table(UMIs), decreasing=TRUE)[[i]]>=2*sort(table(UMIs), decreasing=TRUE)-1)))!=0){
            sortedUMI[[i]] <- rbind(i, which(UMI_levenshtein.out[[i]]<=levenshtein_filt&sort(table(UMIs), decreasing=TRUE)[[i]]>=2*sort(table(UMIs), decreasing=TRUE)-1))
        }else{
            
        }
    }
    
    sortedUMI.bind <- do.call(cbind, sortedUMI)
    UMI_mat_prep <- matrix( sortedUMI.bind, nc = 2, byrow = TRUE)
    
    UMI_graph <- make_graph(UMI_mat_prep, n = length(sortedUMI), directed = FALSE)
    
    return(UMI_graph)
}


.UMIlevenshtein <- function(UMI_unique_sort, levenshtein_max = 6)
{
    UMI_levenshtein <- vector("list", length(UMI_unique_sort))
    
    for (i in 1:length(UMI_unique_sort)){
        UMI_levenshtein[[i]] <- stringdist(UMI_unique_sort[i], UMI_unique_sort)
        UMI_levenshtein[[i]][i] <- levenshtein_max
    } 
    
    return(UMI_levenshtein)
}