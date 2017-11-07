clusterGrouping <- function(overlap, splitted_overlap, grange)

#Returns subgroups of reads that were ordered by their cluster id

{
  
    cluster_id <- rep(0, queryLength(overlap))
  
    for (element in names(splitted_overlap)) {
        reads <- splitted_overlap[[element]]
        change.id <- cluster_id[reads]
        change.id <- change.id[change.id!=0]
        cluster_id[cluster_id%in%change.id] <- as.integer(element)
        cluster_id[reads] <- as.integer(element)
    }
    subgroups <- split(grange, cluster_id) #split reads in respective cluster
    
    return(subgroups)
}
