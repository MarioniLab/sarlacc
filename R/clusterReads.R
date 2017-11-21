clusterReads <- function(map_data, subgroups, geneCount = FALSE, cluster_length = 1)

#Returns list of read names grouped by their cluster id 

{
    cluster.id <- lapply(subgroups, FUN="names") # Read names clustered in same way as reads in subgroups.

    
    gene_count <- .geneCount(grange = map_data$grange, cluster.id = cluster.id)
    cluster.seq <- .seqAssign(map_data = map_data$mapping_data, subgroups = subgroups, cluster.id = cluster.id)

    if(geneCount){
        cluster.seq <- cluster.seq[gene_count==1]
    }
    
    # Filter cluster to contain more reads than cluster_length.
    if(!is.null(cluster_length)){
        len.id <- lapply(cluster.seq, "length")
        clustered.filter <- id_seq_filt[len.id > cluster_length]
        
        return(list(clusterReads = clustered.filter, geneCount = gene_count))
    }else{
        return(list(clusterReads = cluster.seq, geneCount = gene_count))
    }

}



.geneCount <- function(grange, cluster.id)
    #Returns vector with length of number of cluster id's containing the respective number of gene id's
{
    
    #Subgroup gr4 by cluster id reads and store them in cluster_reads as 
    #cluster id subgroups containing respective reads and gene id in mcols
    cluster_reads <- vector('list', length(cluster.id))
    for (element in 1:length(cluster.id)){
        reads <- cluster.id[[element]]
        cluster_reads[[element]] <- grange[reads]
    }
    #Check how many different gene id's are in each cluster id
    gene_count <- vector("integer", length(cluster_reads))
    for (i in 1:length(cluster_reads)){
        gene_count[i] <- length(table(unlist(mcols(cluster_reads[[i]]))))
    }
    
    return(gene_count)
}

.seqAssign <- function(map_data, subgroups, cluster.id)
    #Returns list of cluster id groups that contain the corresponding sequences
{
    seqs <- map_data$SEQ
    names(seqs) <- map_data$QNAME
    
    #List of lists of sequence of reads that are assigned to the same cluster_id
    cluster.seq <- vector("list", length(subgroups))
    for (i in 1:length(subgroups)){
        reads <- cluster.id[[i]]
        cluster.seq[[i]] <- seqs[reads]
    }
    
    return(cluster.seq)
}

