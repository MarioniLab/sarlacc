#' @export
#' @importFrom Biostrings alignedPattern alignedSubject
#' @importFrom BiocGenerics start
#' @importFrom methods as
#' @importClassesFrom IRanges IntegerList
homopolymerMatcher <- function(alignments) {
    ref <- alignedSubject(alignments)
    reads <- alignedPattern(alignments)

    info <- .Call(cxx_match_homopolymers, ref, reads)
    positions <- info[[2]]
    obs.len <- info[[3]]

    output <- homopolymerFinder(ref[1])[[1]]
    ref.pos <- start(output)
    
    m <- match(positions, ref.pos)
    m <- factor(m, levels=seq_along(ref.pos))
    by.pos <- split(obs.len, m, drop=FALSE)
    by.pos <- lapply(by.pos, table)

    by.pos <- as(by.pos, "IntegerList")
    mcols(output)$observed <- by.pos
    return(output)
}
