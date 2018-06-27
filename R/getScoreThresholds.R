#' @export
#' @importFrom S4Vectors metadata
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom BiocParallel SerialParam
getScoreThresholds <- function(aligned, error=0.01, BPPARAM=SerialParam())
# Scrambles the input sequence and performs the same thing as adaptorAlign but with a scrambled input. 
# Identifies the score threshold for the adaptors that achieves the specified error rate.
#
# written by Aaron Lun
# created 10 March 2018
{
    adaptor1 <- metadata(aligned$adaptor1)$sequence
    adaptor2 <- metadata(aligned$adaptor2)$sequence
    tolerance <- metadata(aligned)$tolerance
    go <- metadata(aligned)$gapOpening 
    ge <- metadata(aligned)$gapExtension
    ma <- metadata(aligned)$match
    mm <- metadata(aligned)$mismatch

    # Scrambling the start and end of the read sequences.
    reads.out <- .get_front_and_back(reads, tolerance)
    reads.start <- reads.out$front
    reads.end <- reads.out$back
    
    has.quality <- is(reads, "QualityScaledDNAStringSet")
    scrambled.start <- .scramble_input(reads.start, has.quality)
    scrambled.end <- .scramble_input(reads.end, has.quality)

    all.args <- .setup_alignment_args(has.quality, go, ge, ma, mm)
    scrambled.scores <- .get_all_alignments(adaptor1, adaptor2, scrambled.start, scrambled.end, all.args, scoreOnly=TRUE, BPPARAM=BPPARAM)

    is.reverse <- .resolve_strand(scrambled.scores$start, scrambled.scores$end, 
                                  scrambled.scores$rc.start, scrambled.scores$rc.end)$reversed

    scram.score1 <- ifelse(is.reverse, scrambled.scores$rc.start, scrambled.scores$start)
    scram.score2 <- ifelse(is.reverse, scrambled.scores$rc.end, scrambled.scores$end)
    scram.score1 <- sort(scram.score1)
    scram.score2 <- sort(scram.score2)

    # Computing the expected FDR at each score threshold (a bit conservative, 
    # as we can't remove true positives prior to scrambling).
    score1 <- sort(aligned$adaptor1$score)
    fdr1 <- (length(scram.score1) - findInterval(score1, scram.score1))/(length(score1) - seq_along(score1))
    ix1 <- min(which(fdr1 <= error))

    score2 <- sort(aligned$adaptor2$score)
    fdr2 <- (length(scram.score2) - findInterval(score2, scram.score2))/(length(score2) - seq_along(score2))
    ix2 <- min(which(fdr2 <= error))

    return(list(threshold1=score1[ix1],
                threshold2=score2[ix2],
                scores1=list(reads=score1, scrambled=scram.score1),
                scores2=list(reads=score2, scrambled=scram.score2)))
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
