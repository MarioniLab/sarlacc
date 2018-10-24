#' @export
#' @importFrom Biostrings reverseComplement 
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom BiocParallel bpmapply SerialParam bpstart bpstop bpisup
#' @importFrom BiocGenerics width rownames<-
#' @importFrom ShortRead FastqStreamer yield 
adaptorAlign <- function(adaptor1, adaptor2, filepath, tolerance=250, gapOpening=5, gapExtension=1, 
    qual.type=c("phred", "solexa", "illumina"), number=1e5, BPPARAM=SerialParam())
# This function aligns both adaptors to the read sequence with the specified parameters,
# and returns the alignments that best match the sequence (with reverse complementing if necessary).
#    
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 7 November 2017    
{
    adaptor1 <- toupper(as.character(adaptor1))
    adaptor2 <- toupper(as.character(adaptor2))
    qual.type <- match.arg(qual.type)
    qual.class <- .qual2class(qual.type)

    all.args <- list(adaptor1=adaptor1, adaptor2=adaptor2, tolerance=tolerance, 
            subseq1=.setup_subseqs(adaptor1), subseq2=.setup_subseqs(adaptor2),
            gap.opening=gapOpening, gap.extension=gapExtension)

    # Looping across reads.
    fhandle <- FastqStreamer(filepath, n=number)
    on.exit(close(fhandle))
    all.starts <- all.ends <- all.names <- all.widths <- list()
    counter <- 1L

    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    while (length(fq <- yield(fhandle))) {
        reads <- .FASTQ2QSDS(fq, qual.class)
        by.cores <- .parallelize(reads, BPPARAM)
        out <- bpmapply(by.cores, FUN=.align_AA_internal, MoreArgs=all.args, BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)

        cur.starts <- do.call(rbind, lapply(out, "[[", i="START"))
        cur.ends <- do.call(rbind, lapply(out, "[[", i="END"))
        cur.rc.starts <- do.call(rbind, lapply(out, "[[", i="RSTART"))
        cur.rc.ends <- do.call(rbind, lapply(out, "[[", i="REND"))

        # Standardizing the strand.
        strand.out <- .resolve_strand(cur.starts$score, cur.ends$score, cur.rc.starts$score, cur.rc.ends$score)
        is_reverse <- strand.out$reversed
       
        cur.starts[is_reverse,] <- cur.rc.starts[is_reverse,]
        cur.ends[is_reverse,] <- cur.rc.ends[is_reverse,]

        all.names[[counter]] <- names(reads)
        all.widths[[counter]] <- width(reads)
        all.starts[[counter]] <- cur.starts
        all.ends[[counter]] <- cur.ends

        counter <- counter + 1L
    }

    align_start <- do.call(rbind, all.starts)
    align_end <- do.call(rbind, all.ends)

    details <- list(gapOpening=gapOpening, gapExtension=gapExtension)
    metadata(align_start) <- c(list(sequence=adaptor1), details)
    metadata(align_end) <- c(list(sequence=adaptor2), details)

    # Adjusting the reverse coordinates for the read length.
    old.start <- align_end$start
    old.end <- align_end$end

    all.widths <- unlist(all.widths)
    align_end$start <- all.widths - old.start + 1L
    align_end$end <- all.widths - old.end + 1L

    all.names <- unlist(all.names)
    rownames(align_start) <- rownames(align_end) <- all.names
    output <- DataFrame(read.width=all.widths, adaptor1=I(align_start), adaptor2=I(align_end), reversed=is_reverse, row.names=all.names)
    metadata(output) <- list(filepath=filepath, qual.type=qual.type, tolerance=tolerance)
    return(output)
}

#########################################
#########################################

#' @importFrom Biostrings reverseComplement
#' @importFrom BiocGenerics width
#' @importFrom XVector subseq
.get_front_and_back <- function(reads, tolerance) 
# Extracting the front part of the read (on the forward strand) 
# and the end of the read (on the reverse strand).
{
    tolerance <- pmin(tolerance, width(reads))
    reads.start <- subseq(reads, start = 1, width = tolerance)
    reads.end <- subseq(reads, end = width(reads), width = tolerance)
    reads.end <- reverseComplement(reads.end)
    return(list(front=reads.start, back=reads.end))
}

.qual2class <- function(qual.type) {
     paste0(toupper(substr(qual.type, 1, 1)), substr(qual.type, 2, nchar(qual.type)), "Quality")
}

#' @importFrom Biostrings quality
#' @importFrom methods as
#' @importFrom ShortRead FastqStreamer sread id
.FASTQ2QSDS <- function(fq, qual.class) {
    seq <- sread(fq)
    qual <- as(quality(fq), qual.class)
    reads <- QualityScaledDNAStringSet(seq, qual)
    names(reads) <- as.character(id(fq))
    reads
}

.resolve_strand <- function(start.score, end.score, rc.start.score, rc.end.score) 
# This determines whether the read orientation needs to be flipped so that 
# adaptor 1 is at the start and adaptor 2 is at the end (see ?tuneAlignments 
# for an explanation of what the 'final.score' means).
{
    fscore <- pmax(start.score, 0) + pmax(end.score, 0)
    rscore <- pmax(rc.start.score, 0) + pmax(rc.end.score, 0)
    is.reverse <- fscore < rscore
    final.score <- ifelse(is.reverse, rscore, fscore)
    return(list(reversed=is.reverse, scores=final.score))
}

#' @importFrom BiocParallel SerialParam bpnworkers
#' @importFrom S4Vectors head
.parallelize <- function(reads, BPPARAM) 
# Splits reads across cores. Note the as.list() hack,
# as bplapply doesn't acknowledge DNAStringSetList when BPPARAM is a BatchtoolsParam object.
{ 
    n.cores <- bpnworkers(BPPARAM)
    bounds <- seq(from=1L, to=length(reads), length.out=n.cores+1L)
    ids <- findInterval(seq_along(reads), head(bounds, -1))
    as.list(split(reads, ids))
}

.setup_subseqs <- function(adaptor) {
    mpos <- gregexpr("[^ACTG]+", adaptor)[[1]]
    if (length(mpos)==1L && mpos==-1L) {
        list(starts=integer(0), ends=integer(0))
    } else {
        list(starts=as.integer(mpos), ends=as.integer(mpos) + attr(mpos, "match.length") - 1L)
    }
}

#' @importFrom Biostrings quality pattern subject aligned unaligned
#' @importFrom BiocGenerics score start end
#' @importFrom S4Vectors DataFrame 
#' @importFrom XVector subseq compact
#' @importFrom methods new
#' @importClassesFrom S4Vectors DataFrame
.align_and_extract <- function(adaptor, reads, gap.opening, gap.extension, subseq.starts, subseq.ends)
# Align reads and extract alignment strings.
# This requires custom code as the Biostrings getters for the full alignment strings are not efficient.
{
    quals <- quality(reads)
    out <- .Call(cxx_adaptor_align, reads, quals, .create_encoding_vector(quals), 
            gap.opening, gap.extension, 
            adaptor, subseq.starts - 1L, subseq.ends)
   
    output <- DataFrame(score=out[[1]], start=out[[2]], end=out[[3]])
    
    if (length(subseq.starts)) {
        segments <- vector("list", length(out[[4]]))
        for (i in seq_along(segments)) {
            segments[[i]] <- compact(subseq(reads, out[[4]][[i]], out[[5]][[i]]))
        }
        names(segments) <- sprintf("Sub%i", seq_along(segments))
        segments <- DataFrame(segments)
    } else {
        segments <- new("DataFrame", nrows=length(reads))
    }

    output$subseqs <- segments
    output
}

.align_AA_internal <- function(reads, adaptor1, adaptor2, tolerance, subseq1, subseq2, ...) 
# Wrapper function to pass to bplapply, along with the sarlacc namespace.
{
    reads.out <- .get_front_and_back(reads, tolerance)
    reads.start <- reads.out$front
    reads.end <- reads.out$back

    # Performing the alignment of each adaptor to the start/end of the read.
    list(
        START=.align_and_extract(adaptor=adaptor1, reads=reads.start, subseq.starts=subseq1$starts, subseq.ends=subseq1$ends, ...),
        END=.align_and_extract(adaptor=adaptor2, reads=reads.end, subseq.starts=subseq2$starts, subseq.ends=subseq2$ends, ...),
        RSTART=.align_and_extract(adaptor=adaptor1, reads=reads.end, subseq.starts=subseq1$starts, subseq.ends=subseq1$ends, ...),
        REND=.align_and_extract(adaptor=adaptor2, reads=reads.start, subseq.starts=subseq2$starts, subseq.ends=subseq2$ends, ...)
    )
}
