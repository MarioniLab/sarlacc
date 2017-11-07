clusterGeneReadFilter <- function(id_seq, gene_count)

#Returns a list of cluster id sequences filtered for containing only one gene and more than one read per cluster id

{
    id_seq_filt <- id_seq[gene_count==1] #cluster with one particular gene id
    len.id <- lapply(id_seq_filt, length) #number of reads in each cluster id
    id_seq_len.fil <- id_seq_filt[len.id>1] #filter for more than one read per cluster
  
    return(id_seq_len.fil)
}
