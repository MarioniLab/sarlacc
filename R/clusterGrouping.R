clusterGrouping <- function(overlap.reads, map_data)

# Returns subgroups of reads that were ordered by their cluster id.

{
  
    cluster_id <- rep(0,  as.integer(names(overlap.reads)[length(overlap.reads)]))
  
    for (element in names(overlap.reads)) {
        reads <- overlap.reads[[element]]
        change.id <- cluster_id[reads]
        change.id <- change.id[change.id!=0]
        cluster_id[cluster_id%in%change.id] <- as.integer(element)
        cluster_id[reads] <- as.integer(element)
    }
    cluster_id[cluster_id==0] <- NA
    subgroups <- split(map_data$grange, cluster_id) #split reads in respective cluster
    
    return(subgroups)
}
