#' @export
#' @importFrom Biostrings alignedPattern alignedSubject
#' @importFrom methods as
#' @importClassesFrom IRanges IntegerList
#' @importFrom S4Vectors split DataFrame
errorFinder <- function(alignments) {
    ref <- alignedSubject(alignments)
    reads <- alignedPattern(alignments)

    info <- .Call(cxx_find_errors, ref, reads)
    output <- lapply(info[1:5], function(x) { c(x, NA) }) # adding NA to accommodate one-past-end.
    names(output) <- c("to.A", "to.C", "to.G", "to.T", "deletion")
    output <- DataFrame(output)

    # Aggregating insertion information.
    positions <- factor(info[[6]] + 1L, levels=seq_len(nrow(output)))
    by.pos <- split(info[[7]], positions, drop=FALSE)
    by.pos <- lapply(by.pos, table)
    by.pos <- as(by.pos, "IntegerList")

    total.errors <- Reduce("+", output)
    DataFrame(same=length(reads) - total.errors, output, insertion=by.pos)
}
