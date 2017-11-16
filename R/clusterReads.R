clusterReads <- function(map_data, subgroups)

#Returns list of read names grouped by their cluster id 

{
    clustered.id <- lapply(subgroups, FUN="names") #list of cluster id's containing list of reads that were assigned to the particular cluster id
    
    gene_count <- .geneCount(grange = map_data$grange, clustered.id = clustered.id)
    clustered.seq <- .seqAssign(map_data = map_data$mapping_data, subgroups = subgroups, clustered.id = clustered.id)
    clustered.filtered <- .geneReadFilter(clustered.seq = clustered.seq, gene_count = gene_count)
    
    return(list(Cluster.reads = clustered.filtered, geneCount = gene_count, clusterSeq = clustered.seq))
}



.geneCount <- function(grange, clustered.id)
    #Returns vector with length of number of cluster id's containing the respective number of gene id's
{
    
    #Subgroup gr4 by cluster id reads and store them in cluster_reads as 
    #cluster id subgroups containing respective reads and gene id in mcols
    cluster_reads <- vector('list', length(clustered.id))
    for (element in 1:length(clustered.id)){
        reads <- clustered.id[[element]]
        cluster_reads[[element]] <- grange[reads]
    }
    #Check how many different gene id's are in each cluster id
    gene_count <- vector("integer", length(cluster_reads))
    for (i in 1:length(cluster_reads)){
        gene_count[i] <- length(table(mcols(cluster_reads[[i]])))
    }
    
    return(gene_count)
}

.seqAssign <- function(map_data, subgroups, clustered.id)
    #Returns list of cluster id groups that contain the corresponding sequences
{
    seqs <- map_data$SEQ
    names(seqs) <- map_data$QNAME
    
    #List of lists of sequence of reads that are assigned to the same cluster_id
    clustered.seq <- vector("list", length(subgroups))
    for (i in 1:length(subgroups)){
        reads <- clustered.id[[i]]
        clustered.seq[[i]] <- seqs[reads]
    }
    
    return(clustered.seq)
}

.geneReadFilter <- function(clustered.seq, gene_count)
    #Returns a list of cluster id sequences filtered for containing only one gene and more than one read per cluster id
{
    id_seq_filt <- clustered.seq[gene_count==1] #cluster with one particular gene id
    len.id <- lapply(id_seq_filt, "length") #number of reads in each cluster id
    clustered.filtered <- id_seq_filt[len.id>1] #filter for more than one read per cluster
    
    return(clustered.filtered)
}

