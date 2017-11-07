UMIlevenshtein <- function(UMI_unique_sort, levenshtein_max = 6)

#Filtered UMI's are compared by computing the levenshtein distance to each other UMI available in the list.
#The levenshtein distance is used to measure the difference between two strings.
#It is defined as the minimum number of single character changes (substitutions, insertions, deletions) required to transform one string into the other.
#UMIlevenshtein returns a list containing an integer vector of levenshtein distances for each UMI.

{
    UMI_levenshtein <- vector("list", length(UMI_unique_sort))
    
    for (i in 1:length(UMI_unique_sort)){
        UMI_levenshtein[[i]] <- stringdist(UMI_unique_sort[i], UMI_unique_sort)
        UMI_levenshtein[[i]][i] <- levenshtein_max
    } 
    
    return(UMI_levenshtein)
}
