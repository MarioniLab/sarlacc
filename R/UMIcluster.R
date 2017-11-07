UMIcluster <- function(UMI_unique_sort, UMI_graph)

#UMI's are split corresponding to their graph node.

{
    UMI_clusters <- split(UMI_unique_sort, membership(components(UMI_graph)))
    
    return(UMI_clusters)
}
