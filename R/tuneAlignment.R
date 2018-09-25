#' @export
#' @importFrom Biostrings QualityScaledDNAStringSet
#' @importFrom methods is
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors head
tuneAlignment <- function(adaptor1, adaptor2, reads, tolerance=200, number=10000,
    gapOp.range=c(4, 10), gapExt.range=c(1, 5), match.range=c(1, 2), mismatch.range=c(-1, 0),
    BPPARAM=SerialParam()) 
# This function scrambles the input subject sequences and tries to identify the
# alignment parameters that minimize the overlap between the true alignment 
# scores (to the unscrambled subject) and those to the scrambled sequences. 
#
# written by Aaron Lun
# created 26 February 2018
{
    has.quality <- is(reads, "QualityScaledDNAStringSet")
    adaptor1 <- .assign_qualities(adaptor1, has.quality)
    adaptor2 <- .assign_qualities(adaptor2, has.quality)
    reads <- .assign_qualities(head(reads, number), has.quality)
    if (has.quality) {
        match.range <- mismatch.range <- c(0L, 0L)
    }
    if (length(reads)==0L) {
        return(list(parameters=list(match=NA_integer_, mismatch=NA_integer_, gapOpening=NA_integer_, gapExtension=NA_integer_),
                    scores=list(reads=numeric(0), scrambled=numeric(0))))
    }

    # Taking subsequences of the reads for pairwise alignment.
    reads.out <- .get_front_and_back(reads, tolerance)
    reads.start <- reads.out$front
    reads.end <- reads.out$back
    scrambled.start <- .scramble_input(reads.start, has.quality)
    scrambled.end <- .scramble_input(reads.end, has.quality)

    # Setting up the ranges to be optimized.
    gapOp.range <- as.integer(cummax(gapOp.range))
    gapExt.range <- as.integer(cummax(gapExt.range))
    match.range <- as.integer(cummax(match.range))
    mismatch.range <- as.integer(cummax(mismatch.range))

    # Performing a grid search to maximize the separation in scores.
    max.score <- 0L
    final.match <- final.mismatch <- final.gapOp <- final.gapExt <- NA
    final.reads <- final.scrambled <- NULL

    for (mm in seq(mismatch.range[1], mismatch.range[2], by=1)) { 
        for (ma in seq(match.range[1], match.range[2], by=1)) { 
            for (go in seq(gapOp.range[1], gapOp.range[2], by=1)) { 
                for (ge in seq(gapExt.range[1], gapExt.range[2], by=1)) { 

                    all.args <- .setup_alignment_args(has.quality, go, ge, ma, mm)
                    all.args$BPPARAM <- BPPARAM
                    all.args$scoreOnly <- TRUE

                    align_start <- do.call(.bplalign, c(list(adaptor=adaptor1, reads=reads.start), all.args))
                    align_end <- do.call(.bplalign, c(list(adaptor=adaptor2, reads=reads.end), all.args))
                    align_revcomp_start <- do.call(.bplalign, c(list(adaptor=adaptor1, reads=reads.end), all.args))
                    align_revcomp_end <- do.call(.bplalign, c(list(adaptor=adaptor2, reads=reads.start), all.args))
                    read.scores <- .resolve_strand(align_start, align_end, align_revcomp_start, align_revcomp_end)$scores

                    sc_align_start <- do.call(.bplalign, c(list(adaptor=adaptor1, reads=scrambled.start), all.args))
                    sc_align_end <- do.call(.bplalign, c(list(adaptor=adaptor2, reads=scrambled.end), all.args))
                    sc_align_revcomp_start <- do.call(.bplalign, c(list(adaptor=adaptor1, reads=scrambled.end), all.args))
                    sc_align_revcomp_end <- do.call(.bplalign, c(list(adaptor=adaptor2, reads=scrambled.start), all.args))
                    scrambled.scores <- .resolve_strand(sc_align_start, sc_align_end, sc_align_revcomp_start, sc_align_revcomp_end)$scores

                    cur.score <- .tied_overlap(read.scores, scrambled.scores)
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

    return(list(parameters=list(gapOpening=final.gapOp, gapExtension=final.gapExt, match=final.match, mismatch=final.mismatch),
                scores=list(reads=final.reads, scrambled=final.scrambled)))
}

.tied_overlap <- function(real, fake) 
# Calculating overlap (need average of left.open=TRUE/FALSE to handle ties).
{
    fake <- sort(fake)
    upper.bound <- findInterval(real, fake)
    lower.bound <- findInterval(real, fake, left.open=TRUE)
    sum((upper.bound + lower.bound)/2)/(length(real)*length(fake))
}
