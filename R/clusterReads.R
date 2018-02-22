clusterReads <- function(reads, kmer=5, ndim=50, nclusters=NULL, rand.seed=1000, ...)
# Clusters reads based on their k-mer distances. In particular, 
# generate a shared-nearest neighbour graph of the various reads.
#
# written by Aaron Lun
# created 22 February 2018    
{
    kmer <- as.integer(kmer)
    ndim <- as.integer(ndim)
    mat <- .Call(cxx_get_kmer_matrix, reads, kmer)
    
    # Running irlba to get a low-dimensional approximation.
    if (!is.null(rand.seed)) { 
        set.seed(rand.seed)
    }
    out <- prcomp_irlba(mat, n=ndim)

    # Using FlowSOM to cluster the reads into categories. 
    if (is.null(nclusters)) {
        nclusters <- sqrt(length(reads))
    }
    
    ff <- flowFrame(exprs = out$x)
    fSOM <- ReadInput(ff, compensate = FALSE, transform = FALSE, scale = FALSE, silent = TRUE)
    grid.dim <- ceiling(sqrt(nclusters))
    fSOM <- BuildSOM(fSOM, silent = TRUE, xdim = grid.dim, ydim = grid.dim)
    return(fSOM$map$mapping[,1])
}

