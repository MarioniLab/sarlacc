#' @export
clusterReads <- function(reads, groups, kmer=5, ...)
# Clusters reads based on their k-mer distances. In particular, 
# generate a shared-nearest neighbour graph of the various reads.
#
# written by Aaron Lun
# created 22 February 2018    
{
    kmer <- as.integer(kmer)
    ndim <- as.integer(ndim)
    
    by.group <- split(seq_along(groups), groups)
    output <- integer(length(reads))
    last <- 0L

    for (g in by.group) {
        mat <- .Call(cxx_get_kmer_matrix, reads[g], kmer)
        graph <- scran::buildSNNGraph(mat, transposed=TRUE, ...) 
        P <- igraph::cluster_fast_greedy(graph)
        output[g] <- P + last
        last <- last + 1L
    }

    return(output)
}

