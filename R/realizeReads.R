#' @export
#' @importFrom S4Vectors metadata
#' @importFrom ShortRead FastqStreamer yield
#' @importFrom Biostrings reverseComplement subseq
realizeReads <- function(aligned, number=1e5, trim=TRUE) 
# Converts conceptual reads into actual reads.
#
# written by Aaron Lun
# created 1 October 2018
{
    filepath <- metadata(aligned)$filepath
    qual.type <- metadata(aligned)$qual.type
    qual.class <- .qual2class(qual.type)

    fhandle <- FastqStreamer(filepath, n=number)
    on.exit(close(fhandle))
    all.reads <- list()
    counter <- 1L

    while (length(fq <- yield(fhandle))) {
        reads <- .FASTQ2QSDS(fq, qual.class)
        reads <- reads[names(reads) %in% rownames(aligned)]
        all.reads[[counter]] <- reads
        counter <- counter+1L
    }

    all.reads <- do.call(c, all.reads)
    m <- match(rownames(aligned), names(all.reads))
    if (any(is.na(m))) {
        stop("read names in 'aligned' not present in FASTQ file")
    }
    all.reads <- all.reads[m]

    # Applying reverse-complementing and trimming.
    all.reads[aligned$reversed] <- reverseComplement(all.reads[aligned$reversed])

    if (trim) {
        if (!is.null(aligned$trim.start)) {
            all.reads <- subseq(all.reads, aligned$trim.start, aligned$trim.end)
        } else { 
            warning("no 'trim.start' detected, run 'filterReads' first")
        }
    }

    all.reads
}
