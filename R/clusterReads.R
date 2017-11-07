clusterReads <- function(subgroups)

#Returns list of read names grouped by their cluster id 

{
    id_reads <- lapply(subgroups, names) #list of cluster id's containing list of reads that were assigned to the particular cluster id
    
    return(id_reads)
}
