geneClusterReadAssign <- function(gr_gene, id_reads)

#Returns vector with length of number of cluter id's containing the respective number of gene id's

{
  
    #Subgroup gr4 by cluster id reads and store them in cluster_reads as 
    #cluster id subgroups containing respective reads and gene id in mcols
    cluster_reads <- vector('list', length(id_reads))
    for (element in 1:length(id_reads)){
        reads <- id_reads[[element]]
        cluster_reads[[element]] <-  gr_gene[reads]
    }
  
    #Check how many different gene id's are in each cluster id
    gene_count <- vector("integer", length(cluster_reads))
    for (i in 1:length(cluster_reads)){
        gene_count[i] <- length(table(mcols(cluster_reads[[i]])))
    }
    
    return(gene_count)
}
