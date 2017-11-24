# For use after readOverlap.
# Filters reads for each subgroup.

clusterReads <- function(map_data, overlap.reads, GroupSize = 1)
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
    
    cluster.name <- split(map_data, cluster_id)
    
    greater <- which(unlist(lapply(cluster.name, "length")) >= GroupSize)
    
    cluster.filt <- vector("list", length(greater))
    id.filt <- vector("list", length(greater))
    
    for(i in 1:length(greater)){
        cluster.filt[[i]] <- cluster.name[[greater[i]]]
        id.filt[[i]] <- unique(names(unlist(cluster.name[[greater[i]]])))
    }
    
    # Preparation for msa and consensus sequence.
    # Probably better after UMIgroup.
    # Subgrouping reads with id.
    
    # id.nfilt <- as.integer(id.filt[[1]])
    # id.seq <- chop_data$reads[id.nfilt]
}
    

    
    
    
    
    
  