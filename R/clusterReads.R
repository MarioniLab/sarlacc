#' @export
#' @importFrom Biostrings pairwiseAlignment pid QualityScaledDNAStringSet
#' @importFrom methods is
#' @importFrom igraph make_graph add_edges membership components subcomponent
clusterReads <- function(reads, groups, min.identity=80, gapOpening=5, gapExtension=2, match=1, mismatch=0)
# Clusters reads together based on whether they share sufficient identity.
#
# written by Aaron Lun
# created 22 February 2018    
{
    has.quality <- is(reads, "QualityScaledDNAStringSet")
    align.args <- .setup_alignment_args(has.quality, gapOpening, gapExtension, match, mismatch)
    align.args$type <- "global" # switching to global alignment for comparing reads to each other, especially if they _should_ be the same.

    by.group <- split(seq_along(groups), groups)
    output <- integer(length(reads))
    last <- 0L

    for (gr in by.group) {
        current <- reads[gr]
        all.nodes <- seq_along(gr)
        graph <- make_graph(edges=character(0), isolates=all.nodes, directed=FALSE)

        # Only bother aligning to survivors that aren't already connected.
        for (counter in all.nodes) {
            unconnected <- setdiff(all.nodes, subcomponent(graph, counter))
            aligned <- do.call(pairwiseAlignment, c(list(current[unconnected], current[counter]), align.args))
            keep <- pid(aligned, "PID3") >= min.identity
            if (!any(keep)) {
                next
            }
            graph <- add_edges(graph, rbind(counter, unconnected[keep]))
        }
        
        P <- membership(components(graph))
        output[gr] <- P + last
        last <- last + max(P)
    }

    return(output)
}
