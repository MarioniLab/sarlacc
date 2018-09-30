#' @export
#' @importFrom S4Vectors metadata
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom BiocParallel SerialParam
#' @importFrom ShortRead FastqStreamer yield
getAdaptorThresholds <- function(aligned, error=0.01, block.size=1e8, BPPARAM=SerialParam())
# Scrambles the input sequence and performs the same thing as adaptorAlign but with a scrambled input. 
# Identifies the score threshold for the adaptors that achieves the specified error rate.
#
# written by Aaron Lun
# created 10 March 2018
{
    adaptor1 <- metadata(aligned$adaptor1)$sequence
    adaptor2 <- metadata(aligned$adaptor2)$sequence

    go <- metadata(aligned$adaptor1)$gapOpening 
    ge <- metadata(aligned$adaptor1)$gapExtension
    all.args <- .setup_alignment_args(TRUE, go, ge)
    all.args$BPPARAM <- BPPARAM
    all.args$scoreOnly <- TRUE

    tolerance <- metadata(aligned)$tolerance
    filepath <- metadata(aligned)$filepath
    qual.type <- metadata(aligned)$qual.type
    qual.class <- .qual2class(qual.type)

    # Looping through the files.
    fhandle <- FastqStreamer(filepath, readerBlockSize=block.size)
    on.exit(close(fhandle))
    scram.score1 <- scram.score2 <- used.names <- list()
    counter <- 1L

    while (length(fq <- yield(fhandle))) {
        reads <- .FASTQ2QSDS(fq, qual.class)
        reads <- reads[names(reads) %in% rownames(aligned)]
        used.names[[counter]] <- names(reads)

        # Scrambling the start and end of the read sequences.
        reads.out <- .get_front_and_back(reads, tolerance)
        reads.start <- reads.out$front
        reads.end <- reads.out$back

        scrambled.start <- .scramble_input(reads.start, TRUE)
        scrambled.end <- .scramble_input(reads.end, TRUE)

        # Computing alignment scores to the adaptors.
        scrambled_start <- do.call(.bplalign, c(list(adaptor=adaptor1, reads=scrambled.start), all.args))
        scrambled_end <- do.call(.bplalign, c(list(adaptor=adaptor2, reads=scrambled.end), all.args))
        scrambled_revcomp_start <- do.call(.bplalign, c(list(adaptor=adaptor1, reads=scrambled.end), all.args))
        scrambled_revcomp_end <- do.call(.bplalign, c(list(adaptor=adaptor2, reads=scrambled.start), all.args))

        is.reverse <- .resolve_strand(scrambled_start, scrambled_end, scrambled_revcomp_start, scrambled_revcomp_end)$reversed
        scram.score1[[counter]] <- ifelse(is.reverse, scrambled_revcomp_start, scrambled_start)
        scram.score2[[counter]] <- ifelse(is.reverse, scrambled_revcomp_end, scrambled_end)
        counter <- counter + 1L
    }

    scram.score1 <- unlist(scram.score1)
    scram.score2 <- unlist(scram.score2)
    m <- match(unlist(used.names), rownames(aligned))
    real.score1 <- aligned$adaptor1$score[m]
    real.score2 <- aligned$adaptor2$score[m]

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
