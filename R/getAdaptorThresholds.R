#' @export
#' @importFrom S4Vectors metadata
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom BiocParallel SerialParam
getAdaptorThresholds <- function(aligned, error=0.01, BPPARAM=SerialParam())
# Scrambles the input sequence and performs the same thing as adaptorAlign but with a scrambled input. 
# Identifies the score threshold for the adaptors that achieves the specified error rate.
#
# written by Aaron Lun
# created 10 March 2018
{
    reads <- aligned$reads
    adaptor1 <- metadata(aligned$adaptor1)$sequence
    adaptor2 <- metadata(aligned$adaptor2)$sequence
    tolerance <- metadata(aligned$adaptor1)$tolerance
    go <- metadata(aligned$adaptor1)$gapOpening 
    ge <- metadata(aligned$adaptor1)$gapExtension
    ma <- metadata(aligned$adaptor1)$match
    mm <- metadata(aligned$adaptor1)$mismatch

    # Scrambling the start and end of the read sequences.
    reads.out <- .get_front_and_back(reads, tolerance)
    reads.start <- reads.out$front
    reads.end <- reads.out$back
    
    has.quality <- is(reads, "QualityScaledDNAStringSet")
    scrambled.start <- .scramble_input(reads.start, has.quality)
    scrambled.end <- .scramble_input(reads.end, has.quality)

    # Computing alignment scores to the adaptors.
    all.args <- .setup_alignment_args(has.quality, go, ge, ma, mm)
    all.args$BPPARAM <- BPPARAM
    all.args$scoreOnly <- TRUE

    scrambled_start <- do.call(.bplalign, c(list(adaptor=adaptor1, reads=scrambled.start), all.args))
    scrambled_end <- do.call(.bplalign, c(list(adaptor=adaptor2, reads=scrambled.end), all.args))
    scrambled_revcomp_start <- do.call(.bplalign, c(list(adaptor=adaptor1, reads=scrambled.end), all.args))
    scrambled_revcomp_end <- do.call(.bplalign, c(list(adaptor=adaptor2, reads=scrambled.start), all.args))
    is.reverse <- .resolve_strand(scrambled_start, scrambled_end, scrambled_revcomp_start, scrambled_revcomp_end)$reversed

    scram.score1 <- ifelse(is.reverse, scrambled_revcomp_start, scrambled_start)
    scram.score2 <- ifelse(is.reverse, scrambled_revcomp_end, scrambled_end)
    real.score1 <- aligned$adaptor1$score
    real.score2 <- aligned$adaptor2$score

    list(
        threshold1=.compute_threshold(real.score1, scram.score1, error),
        threshold2=.compute_threshold(real.score2, scram.score2, error),
        scores1=list(reads=real.score1, scrambled=scram.score1),
        scores2=list(reads=real.score2, scrambled=scram.score2)
    )
}

#' @importFrom Biostrings DNAStringSet QualityScaledDNAStringSet
.scramble_input <- function(seqs, has.qual) 
# Scrambles the input sequences. Uses R loops,
# but it should be fast enough for our purposes.   
{
    collected.seqs <- strsplit(as.character(seqs), "")
    if (has.qual) { 
        collected.quals <- strsplit(as.character(quality(seqs)), "")
    }

    for (i in seq_along(seqs)) {
        current <- collected.seqs[[i]]
        o <- sample(length(current))
        collected.seqs[[i]] <- paste(current[o], collapse="")

        if (has.qual) { 
            collected.quals[[i]] <- paste(collected.quals[[i]][o], collapse="")
        }
    }

    output <- DNAStringSet(unlist(collected.seqs))
    if (has.qual) {
        output <- QualityScaledDNAStringSet(output, as(unlist(collected.quals), class(quality(seqs))))
    }
    return(output)
}

.compute_threshold <- function(real, scrambled, error) 
# Computing the expected FDR at each score threshold 
# (a bit conservative, s we can't remove true positives prior to scrambling)
# and then finds the score threshold where the FDR is below the specified error threshold.
{
    real <- sort(real)
    scrambled <- sort(scrambled)
    fdr <- (length(scrambled) - findInterval(real, scrambled))/(length(real) - seq_along(real))
    real[min(which(fdr <= error))]
}
