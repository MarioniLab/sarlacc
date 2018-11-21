#' @export
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom igraph make_graph cluster_fast_greedy
#' @importFrom BiocGenerics width intersect
groupAligned <- function(locations, threshold=0.8, cluster.fun=cluster_fast_greedy) 
# Groups alignments based on their position and overlap.
#
# written by Aaron Lun
# created 22 November 2018
{
    olap <- findOverlaps(locations, ignore.strand=FALSE)
    all.q <- queryHits(olap)
    all.s <- subjectHits(olap)

    # Remove redundant pairs in self-overlap.
    keep <- all.q < all.s
    all.q <- all.q[keep]
    all.s <- all.s[keep]

    # Checking the span of the overlaps.
    left <- locations[all.q]
    right <- locations[all.s]
    overlapping <- intersect(left, right)

    olap.width <- vapply(width(overlapping), sum, FUN.VALUE=0L)
    left.width <- vapply(width(left), sum, FUN.VALUE=0L)
    right.width <- vapply(width(right), sum, FUN.VALUE=0L)
    olap.prop <- olap.width / pmax(left.width, right.width)

    # Removing poor overlaps.
    keep <- olap.prop >= threshold
    edges <- rbind(all.q[keep], all.s[keep])
    g <- make_graph(edges, n=length(locations), directed=FALSE)
    
    out <- cluster.fun(g)$membership
    names(out) <- names(locations)
    out
}
