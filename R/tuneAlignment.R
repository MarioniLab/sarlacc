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
    pre.out <- .preprocess_input(adaptor1, adaptor2, reads)
    adaptor1 <- pre.out$adaptor1
    adaptor2 <- pre.out$adaptor2
    reads <- pre.out$reads

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

                    all.read.scores <- .get_all_alignments(adaptor1, adaptor2, reads.start, reads.end, all.args, scoreOnly=TRUE)
                    read.scores <- .resolve_strand(all.read.scores$start, all.read.scores$end, 
                                                   all.read.scores$rc.start, all.read.scores$rc.end)$scores

                    all.scrambled.scores <- .get_all_alignments(adaptor1, adaptor2, scrambled.start, scrambled.end, all.args, scoreOnly=TRUE)
                    scrambled.scores <- .resolve_strand(all.scrambled.scores$start, all.scrambled.scores$end, 
                                                        all.scrambled.scores$rc.start, all.scrambled.scores$rc.end)$scores
                    
                    cur.score <- sum(findInterval(read.scores, sort(scrambled.scores)))
                    if (max.score < cur.score) {
                        max.score <- cur.score

                        final.gapOp <- go
                        final.gapExt <- ge
                        if (!has.quality) {
                            final.match <- ma
                            final.mismatch <- mm
                        }

                        final.reads <- read.scores
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
