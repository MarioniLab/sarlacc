UMIfilt <- function(UMIs)

#removes duplciated UMI's and sorts them by occurence

{
    if(length(UMIs)>1){
        UMI1_unique_sort <- DNAStringSet(names(sort(table(UMIs$UMI1), decreasing=TRUE))) 
        UMI2_unique_sort <- DNAStringSet(names(sort(table(UMIs$UMI2), decreasing=TRUE)))
        
        return(list(UMI1unique = UMI1_unique_sort, UMI2unique = UMI2_unique_sort))
    } else{
        UMI_unique_sort <- DNAStringSet(names(sort(table(UMIs[[1]]), decreasing=TRUE))) 
        
        return(UMIunique = UMI_unique_sort)
    }
}

