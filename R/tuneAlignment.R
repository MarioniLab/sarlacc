#' @export
#' @importFrom Biostrings QualityScaledDNAStringSet
#' @importFrom methods is
tuneAlignment <- function(adaptor1, adaptor2, reads, tolerance=100, 
                          gapOp.range=c(4, 10), gapExt.range=c(1, 5), 
                          match.range=c(1, 2), mismatch.range=c(-1, 0)) 
# This function scrambles the input subject sequences and tries to identify the
# alignment parameters that minimize the overlap between the true alignment 
# scores (to the unscrambled subject) and those to the scrambled sequences. 
#
# written by Aaron Lun
# created 26 February 2018
{
    if (is.character(adaptor1)) {
        adaptor1 <- DNAString(adaptor1)
    } 
    if (is.character(adaptor2)) { 
        adaptor2 <- DNAString(adaptor2)
    }
    if (is.character(reads)) {
        reads <- DNAStringSet(reads)
    }

    # Taking subsequences of the reads for pairwise alignment.
    reads.out <- .get_front_and_back(reads, tolerance)
    reads.start <- reads.out$front
    reads.end <- reads.out$back

    scrambled.start <- .scramble_input(reads.start)
    scrambled.end <- .scramble_input(reads.end)

    # Setting up the ranges to be optimized.
    gapOp.range <- as.integer(cummax(gapOp.range))
    gapExt.range <- as.integer(cummax(gapExt.range))
    match.range <- as.integer(cummax(match.range))
    mismatch.range <- as.integer(cummax(mismatch.range))

    has.quality <- is(reads, "QualityScaledDNAStringSet")
    if (has.quality) {
        match.range <- mismatch.range <- c(0L, 0L)
    }

    # Performing a grid search to maximize the separation in scores.
    max.score <- 0L
    final.match <- final.mismatch <- final.gapOp <- final.gapExt <- NA
    final.reads <- final.scrambled <- NULL

    for (mm in seq(mismatch.range[1], mismatch.range[2], by=1)) { 
        for (ma in seq(match.range[1], match.range[2], by=1)) { 
            for (go in seq(gapOp.range[1], gapOp.range[2], by=1)) { 
                for (ge in seq(gapExt.range[1], gapExt.range[2], by=1)) { 
        
                    all.args <- .setup_alignment_args(has.quality, go, ge, ma, mm)
                    reads.scores <- .obtain_max_scores(adaptor1, adaptor2, reads.start, reads.end, all.args)
                    scrambled.scores <- .obtain_max_scores(adaptor1, adaptor2, scrambled.start, scrambled.end, all.args)
                    
                    cur.score <- sum(findInterval(reads.scores, sort(scrambled.scores)))
                    if (max.score < cur.score) {
                        max.score <- cur.score

                        final.gapOp <- go
                        final.gapExt <- ge
                        if (!has.quality) {
                            final.match <- ma
                            final.mismatch <- mm
                        }

                        final.reads <- reads.scores
                        final.scrambled <- scrambled.scores
                    }
                }
            }
        }
    }

    return(list(parameters=list(gapOpening=final.gapOp, gapExtension=final.gapExt, 
                                match=final.match, mismatch=final.mismatch),
                scores=list(reads=final.reads, scrambled=final.scrambled)))
}

.scramble_input <- function(seqs) 
# Scrambles the input sequences. Uses R loops,
# but it should be fast enough for our purposes.   
{
    for (i in seq_along(seqs)) {
        current <- seqs[[i]]
        seqs[[i]] <- current[sample(length(current))]
    }
    return(seqs)
}

#' @importFrom Biostrings pairwiseAlignment
.obtain_max_scores <- function(adaptor1, adaptor2, reads.start, reads.end, all.args=list()) {
    align_start <- do.call(pairwiseAlignment, c(list(pattern=reads.start, subject=adaptor1, scoreOnly=TRUE), all.args))
    align_end <- do.call(pairwiseAlignment, c(list(pattern=reads.end, subject=adaptor2, scoreOnly=TRUE), all.args))
    align_revcomp_start <- do.call(pairwiseAlignment, c(list(pattern=reads.end, subject=adaptor1, scoreOnly=TRUE), all.args))
    align_revcomp_end <- do.call(pairwiseAlignment, c(list(pattern=reads.start, subject=adaptor2, scoreOnly=TRUE), all.args))
    pmax(align_start, align_end, align_revcomp_start, align_revcomp_end)
}
