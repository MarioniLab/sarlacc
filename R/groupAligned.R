#' @export
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom igraph make_graph cluster_fast_greedy
#' @importFrom BiocGenerics width intersect
#' @importFrom methods is
#' @importClassesFrom GenomicRanges GenomicRanges
groupAligned <- function(locations, threshold=0.8, cluster.fun=cluster_fast_greedy) 
# Groups alignments based on their position and overlap.
#
# written by Aaron Lun
# created 22 November 2018
{
    mapped <- seqnames(locations)!="*"
    olap <- findOverlaps(locations[mapped], ignore.strand=FALSE)
    all.q <- queryHits(olap)
    all.s <- subjectHits(olap)

    # Adjust for changing of indices when subsetting a GR (but not GRL).
    if (is(locations, "GenomicRanges")) {
        mapped <- which(mapped)
        all.q <- mapped[all.q]
        all.s <- mapped[all.s]
    }

    # Remove redundant pairs, unmapped reads in self-overlap.
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

    keep <- olap.prop >= threshold
    all.q <- all.q[keep]
    all.s <- all.s[keep]

    # Creating the graph and clustering.
    edges <- rbind(all.q, all.s)
    g <- make_graph(edges, n=length(locations), directed=FALSE)
    out <- cluster.fun(g)$membership

    if (is(locations, "GenomicRanges")) {
        out[-mapped] <- NA_integer_
    } else {
        out[!any(mapped)] <- NA_integer_
    }

    names(out) <- names(locations)
    out
}
