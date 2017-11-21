clusterReads <- function(map_data, overlap.reads, geneCount = FALSE, cluster_length = 1)

#Returns list of read names grouped by their cluster id 

{
    id <- .clusterGrouping(overlap.reads = overlap.reads, map_data = map_data)
    subgroups <- split(map_data$grange, id) #split reads in respective cluster
    
    cluster.name <- lapply(subgroups, FUN="names") # Read names clustered in same way as reads in subgroups.
    
    gene_count <- .geneCount(grange = map_data$grange, cluster.name = cluster.name)
    cluster.seq <- .seqAssign(map_data = map_data$mapping_data, subgroups = subgroups, cluster.name = cluster.name)

    if(geneCount){
        cluster.seq <- cluster.seq[gene_count==1]
    }
    
    # Filter cluster to contain more reads than cluster_length.
    if(!is.null(cluster_length)){
        len.id <- lapply(cluster.seq, "length")
        clustered.filter <- cluster.seq[len.id > cluster_length]
    }
    
    if(!is.null(cluster_length)){
        return(list(clusterReads = clustered.filter, geneCount = gene_count, subgroups = subgroups))
    }else{
        return(list(clusterReads = cluster.seq, geneCount = gene_count, subgroups = subgroups))
    }

}



.clusterGrouping <- function(overlap.reads, map_data)
    
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
    
    
    return(id = cluster_id)
}




.geneCount <- function(grange, cluster.name)
    #Returns vector with length of number of cluster id's containing the respective number of gene id's
{
    
    #Subgroup gr4 by cluster id reads and store them in cluster_reads as 
    #cluster id subgroups containing respective reads and gene id in mcols
    cluster_reads <- vector('list', length(cluster.name))
    for (element in 1:length(cluster.name)){
        reads <- cluster.name[[element]]
        cluster_reads[[element]] <- grange[reads]
    }
    #Check how many different gene id's are in each cluster id
    gene_count <- vector("integer", length(cluster_reads))
    for (i in 1:length(cluster_reads)){
        gene_count[i] <- length(table(unlist(mcols(cluster_reads[[i]]))))
    }
    
    return(gene_count)
}

.seqAssign <- function(map_data, subgroups, cluster.name)
    #Returns list of cluster id groups that contain the corresponding sequences
{
    seqs <- map_data$SEQ
    names(seqs) <- map_data$QNAME
    
    #List of lists of sequence of reads that are assigned to the same cluster_id
    cluster.seq <- vector("list", length(subgroups))
    for (i in 1:length(subgroups)){
        reads <- cluster.name[[i]]
        cluster.seq[[i]] <- seqs[reads]
    }
    
    return(cluster.seq)
}

