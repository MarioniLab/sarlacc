UMIindex <- function(UMI_clusters, UMIs)

#Initial UMI's (before sorting and filtering for duplicates) are indexed to receive an useful index for the corresponding reads (sequences).
#UMIindex takes the initial corresponding UMI-read index into account and assigns each UMI node the index numbers of the initial UMI's that are part of the node by either being a duplicate or being connected to the id by their alignment.
#This final list of index numbers makes it possible to subgroup the reads by their UMI's.

{
    UMI_index <- vector("list", length(UMI_clusters))
    for (i in 1:length(UMI_clusters)){
        len <- 1
        for (a in 1:length(UMI_clusters[[i]])){
            id <- which(UMIs==UMI_clusters[[i]][a])
            for (b in 1:length(id)){
                UMI_index[[i]][len] <- id[b] 
                len <- len+1
            }
        }
    }
    
    return(UMI_index)
}
