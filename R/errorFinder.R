#' @export
#' @importFrom Biostrings alignedPattern alignedSubject
#' @importFrom methods as is
#' @importClassesFrom IRanges RleList
#' @importClassesFrom Biostrings DNAStringSet
#' @importFrom S4Vectors split DataFrame Rle
errorFinder <- function(alignments) {
    ref <- alignedSubject(alignments)
    reads <- alignedPattern(alignments)
    if (!is(ref, "DNAStringSet") || !is(reads, "DNAStringSet")) {
        stop("alignments should involve DNAString(Set) objects");
    }

    info <- .Call(cxx_find_errors, ref, reads)
    output <- lapply(info[1:6], function(x) { c(x, NA) }) # adding NA to accommodate one-past-end.
    names(output) <- c("base", "A", "C", "G", "T", "deletion")
    output <- DataFrame(output)

    # Aggregating insertion information.
    positions <- factor(info[[7]] + 1L, levels=seq_len(nrow(output)))
    by.pos <- split(info[[8]], positions, drop=FALSE)
    by.pos <- lapply(by.pos, FUN=function(ilen) {
        tab <- table(ilen)
        Rle(
            c(0L, as.integer(names(tab))), # insertion lengths (adding zeroes for completeness).
            c(length(alignments) - length(ilen), tab) # frequencies of insertion lengths.
        )
    })
    names(by.pos) <- NULL
    by.pos <- as(by.pos, "RleList")

    # Compiling relevant statistics.
    full.stats <- DataFrame(output, insertion=by.pos)
    transitions <- matrix(0L, 4, 4)
    rownames(transitions) <- colnames(transitions) <- c("A", "C", "G", "T") 
    for (base in rownames(transitions)) {
        current <- full.stats[which(full.stats$base==base),colnames(transitions)]
        transitions[base,] <- vapply(current, sum, FUN.VALUE=0L)
    }
    return(list(full=full.stats, transition=transitions))
}
