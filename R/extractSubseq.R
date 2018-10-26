#' @export
#' @importFrom S4Vectors mcols<- metadata
#' @importFrom ShortRead FastqStreamer yield
#' @importFrom BiocParallel bpmapply SerialParam bpstart bpstop bpisup
extractSubseq <- function(aligned, subseq1, subseq2, number=1e5, BPPARAM=SerialParam()) 
# Extracts arbitrary subsequences from the first and second adaptors.
# This requires a realignment as the full alignment info is not held.
#
# written by Aaron Lun
# created 26 October 2018
{
    filepath <- metadata(aligned)$filepath
    qual.type <- metadata(aligned)$qual.type

    go <- metadata(aligned$adaptor1)$gapOpening 
    ge <- metadata(aligned$adaptor1)$gapExtension

    adaptor1 <- metadata(aligned$adaptor1)$sequence
    adaptor2 <- metadata(aligned$adaptor2)$sequence
    tolerance <- metadata(aligned)$tolerance

    # At least one of these should be specified.
    do1 <- do2 <- TRUE
    if (missing(subseq2)) {
        do2 <- FALSE
        subseq2 <- list(starts=integer(0), ends=integer(0))
    } else if (missing(subseq1)) {
        do1 <- FALSE
        subseq1 <- list(starts=integer(0), ends=integer(0))
    }

    qual.class <- .qual2class(qual.type)
    fhandle <- FastqStreamer(filepath, n=number)
    on.exit(close(fhandle))
    all.adaptor1 <- all.adaptor2 <- list()
    counter <- 1L
   
    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    while (length(fq <- yield(fhandle))) {
        reads <- .FASTQ2QSDS(fq, qual.class)

        m <- match(names(reads), rownames(aligned))
        keep <- !is.na(m)
        reads <- reads[keep]
        m <- m[keep]
        mcols(reads)$reversed <- aligned$reversed[m]

        by.cores <- .parallelize(reads, BPPARAM)
        out <- bpmapply(by.cores, FUN=.extract_internal, 
            MoreArgs=list(adaptor1=adaptor1, adaptor2=adaptor2, tolerance=tolerance, 
                subseq1=subseq1, subseq2=subseq2, gap.opening=go, gap.extension=ge),
            BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)

        # Sanity check on scores to ensure the alignments are the same as reported.
        if (do1) {
            all.1.score <- unlist(lapply(out, FUN="[[", i="adaptor1.score"))
            if (!isTRUE(all.equal(all.1.score, aligned$adaptor1$score[m]))) {
                stop("score mismatch from 'aligned' for adaptor 1")
            }
            all.adaptor1[[counter]] <- lapply(out, "[[", i="adaptor1.subseq")
        }

        if (do2) {
            all.2.score <- unlist(lapply(out, FUN="[[", i="adaptor2.score"))
            if (!isTRUE(all.equal(all.2.score, aligned$adaptor2$score[m]))) {
                stop("score mismatch from 'aligned' for adaptor 2")
            }
            all.adaptor2[[counter]] <- lapply(out, "[[", i="adaptor2.subseq")
        }

        counter <- counter+1L
    }

    output <- list()
    if (do1) {
        output$adaptor1 <- do.call(rbind, unlist(all.adaptor1))
    } 
    if (do2) {
        output$adaptor2 <- do.call(rbind, unlist(all.adaptor2))
    }
    output
}

#' @importFrom S4Vectors mcols mcols<-
.extract_internal <- function(reads, adaptor1, adaptor2, tolerance, subseq1, subseq2, ...) 
# Wrapper function to pass to bplapply, along with the sarlacc namespace.
{
    flipped <- mcols(reads)$reversed
    mcols(reads) <- NULL
    reads.out <- .get_front_and_back(reads, tolerance)
    reads.start <- reads.out$front
    reads.end <- reads.out$back

    # Reduce the number of alignments by using existing orientation info.
    actual.starts <- reads.start
    actual.starts[flipped] <- reads.end[flipped]
    actual.ends <- reads.end
    actual.ends[flipped] <- reads.start[flipped]

    if (any(lengths(subseq1)!=0L)) {
        adaptor1 <- .align_and_extract(adaptor=adaptor1, reads=actual.starts, subseq.starts=subseq1$starts, subseq.ends=subseq1$ends, ...)
    } else {
        adaptor1 <- NULL
    }
    if (any(lengths(subseq2)!=0L)) {
        adaptor2 <- .align_and_extract(adaptor=adaptor2, reads=actual.ends, subseq.starts=subseq2$starts, subseq.ends=subseq2$ends, ...)
    } else {
        adaptor2 <- NULL
    }

    return(list(adaptor1.score=adaptor1$score, adaptor1.subseq=adaptor1$subseq, 
                adaptor2.score=adaptor2$score, adaptor2.subseq=adaptor2$subseq)) 
}

