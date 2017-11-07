UMIfilt <- function(UMIs)

#removes duplciated UMI's and sorts them by occurence

{
    UMI_unique_sort <- DNAStringSet(names(sort(table(UMIs), decreasing=TRUE))) 
    
    return(UMI_unique_sort)
}
