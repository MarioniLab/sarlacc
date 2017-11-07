UMImatrix <- function(UMI_levenshtein, UMI_unique_sort, UMIs, levenshtein_filt = 1)

#UMImatrix takes the levenshtein distances into account.
#It sets 1 for UMI's that have a levenshtein distance <=1 and match the equation na >= 2*nb-1. 
#Whereas na is the first/highest count node available and nb is the second highest count node.

{
    UMImat <- matrix(0, length(UMI_levenshtein), length(UMI_levenshtein))
    rownames(UMImat) <- colnames(UMImat) <- as.character(UMI_unique_sort)

    for(i in 1:(length(UMI_unique_sort)-1)){
        if (length(which(UMI_levenshtein[[i]]<=levenshtein_filt&sort(table(UMIs), decreasing=TRUE)[[i]]>=2*sort(table(UMIs), decreasing=TRUE)-1))!=0){
            len <- length(which(UMI_levenshtein[[i]]<=levenshtein_filt&sort(table(UMIs), decreasing=TRUE)[[i]]>=2*sort(table(UMIs), decreasing=TRUE)-1))
            for (a in 1:len){
                UMImat[i,which(UMI_levenshtein[[i]]<=levenshtein_filt&sort(table(UMIs), decreasing=TRUE)[[i]]>=2*sort(table(UMIs), decreasing=TRUE)-1)[a]] <- 1
            }
        
        }else{
        }
    }
    
    return(UMImat)
}
