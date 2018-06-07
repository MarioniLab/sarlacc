#' @export
#' @importFrom methods is
#' @importClassesFrom Biostrings DNAStringSet
#' @importFrom S4Vectors split mcols
#' @importFrom IRanges IRanges
homopolymerFinder <- function(seq)
# Finds homopolymers in a given DNAStringSet object.
{
    if (!is(seq, "DNAStringSet")) {
        stop("'seq' should be a DNAStringSet object")
    }

    all.homo <- .Call(cxx_find_homopolymers, seq)
    homo.range <- IRanges(all.homo[[2]], all.homo[[2]] + all.homo[[3]] - 1L)
    mcols(homo.range)$base <- all.homo[[4]]

    groupings <- factor(all.homo[[1]]+1L, levels=seq_along(seq))
    output <- split(homo.range, groupings, drop=FALSE)
    names(output) <- names(seq)
    return(output)
}

