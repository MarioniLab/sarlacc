#' @export
#' @importFrom Biostrings reverseComplement 
#' @importFrom S4Vectors DataFrame metadata<-
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom methods is
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocGenerics width rownames<-
#' @importFrom ShortRead FastqStreamer yield 
adaptorAlign <- function(adaptor1, adaptor2, filepath, tolerance=250, gapOpening=5, gapExtension=1, 
    qual.type=c("phred", "solexa", "illumina"), block.size=1e8, BPPARAM=SerialParam())
# This function aligns both adaptors to the read sequence with the specified parameters,
# and returns the alignments that best match the sequence (with reverse complementing if necessary).
#    
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 7 November 2017    
{
    adaptor1 <- .assign_qualities(adaptor1, TRUE)
    adaptor2 <- .assign_qualities(adaptor2, TRUE)
    qual.type <- match.arg(qual.type)
    qual.class <- .qual2class(qual.type)

    all.args <- .setup_alignment_args(TRUE, gapOpening, gapExtension)
    all.args$BPPARAM <- BPPARAM

    # Looping across reads.
    fhandle <- FastqStreamer(filepath, readerBlockSize=block.size)
    on.exit(close(fhandle))

    all.starts <- all.ends <- all.rc.starts <- all.rc.ends <- all.names <- all.widths <- list()
    counter <- 1L

    while (length(fq <- yield(fhandle))) {
        reads <- .FASTQ2QSDS(fq, qual.class)
        all.names[[counter]] <- names(reads)
        all.widths[[counter]] <- width(reads)

        # Getting the start and (rc'd) end of the read.
        reads.out <- .get_front_and_back(reads, tolerance)
        reads.start <- reads.out$front
        reads.end <- reads.out$back
    
        # Performing the alignment of each adaptor to the start/end of the read.
        all.starts[[counter]] <- do.call(.bplalign, c(list(adaptor=adaptor1, reads=reads.start), all.args))
        all.ends[[counter]] <- do.call(.bplalign, c(list(adaptor=adaptor2, reads=reads.end), all.args))
        all.rc.starts[[counter]] <- do.call(.bplalign, c(list(adaptor=adaptor1, reads=reads.end), all.args))
        all.rc.ends[[counter]] <- do.call(.bplalign, c(list(adaptor=adaptor2, reads=reads.start), all.args))
        counter <- counter + 1L
    }

    align_start <- do.call(rbind, all.starts)
    align_end <- do.call(rbind, all.ends)
    align_revcomp_start <- do.call(rbind, all.rc.starts)
    align_revcomp_end <- do.call(rbind, all.rc.ends)

    # Standardizing the strand.
    strand.out <- .resolve_strand(align_start$score, align_end$score, align_revcomp_start$score, align_revcomp_end$score)
    is_reverse <- strand.out$reversed

    align_start[is_reverse,] <- align_revcomp_start[is_reverse,]
    align_end[is_reverse,] <- align_revcomp_end[is_reverse,]

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

############################################
########### Internal functions #############
############################################

#' @importFrom Biostrings DNAStringSet PhredQuality QualityScaledDNAStringSet
#' @importFrom BiocGenerics width
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom methods is
.assign_qualities <- function(truth, add.quality=TRUE)
# pairwiseAlignment() needs both strings to have qualities for quality-based pairwise alignment. 
# Obviously, the true reference has no qualities so we assign it the highest score possible.
{
    has.qual <- is(truth, "QualityScaledDNAStringSet")
    if (is.character(truth)) {
        truth <- DNAStringSet(truth)
    } 
    if (add.quality && !has.qual){
        enc <- as.character(PhredQuality(100L))
        qual <- PhredQuality(strrep(enc,  width(truth)))
        truth <- QualityScaledDNAStringSet(truth, qual)
    }
	truth  
}

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

#' @importFrom Biostrings nucleotideSubstitutionMatrix
.setup_alignment_args <- function(has.quality, gapOpening, gapExtension, match, mismatch, type="local-global") {
    all.args <- list(type=type, gapOpening=gapOpening, gapExtension=gapExtension)
    if (!has.quality) { 
        all.args$substitutionMatrix <- nucleotideSubstitutionMatrix(match=match, mismatch=mismatch)
    } else {
        all.args$fuzzyMatrix <- nucleotideSubstitutionMatrix() # support IUPAC ambiguous codes in quality alignments. 
    }
    return(all.args)
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

#' @importFrom Biostrings pairwiseAlignment quality pattern subject aligned unaligned
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#' @importFrom methods is
#' @importFrom BiocGenerics score start end
#' @importFrom S4Vectors DataFrame head
#' @importFrom XVector subseq
.bplalign <- function(adaptor, reads, ..., scoreOnly=FALSE, BPPARAM=SerialParam())
# Align reads in parallel.
{
    n.cores <- bpnworkers(BPPARAM)
    bounds <- seq(from=1L, to=length(reads), length.out=n.cores+1L)
    ids <- findInterval(seq_along(reads), head(bounds, -1))
    by.core <- split(reads, ids)

    by.core <- as.list(by.core) # hack as bplapply doesn't acknowledge DNAStringSetList when BPPARAM is a BatchtoolsParam object.
    out <- bplapply(by.core, FUN=pairwiseAlignment, subject=adaptor, ..., scoreOnly=scoreOnly, BPPARAM=BPPARAM)
    if (scoreOnly) {
        if (length(out)==0L) {
            return(numeric(0))
        }
        return(unlist(out, use.names=FALSE))
    }

    has.quality <- is(reads, "QualityScaledDNAStringSet")
    collected <- vector("list", n.cores)
    collected[[1]] <- DataFrame(score=numeric(0), adaptor=character(0), read=character(0), start=integer(0), end=integer(0)) # Defined output even if there are no reads.
    if (has.quality) {
        collected[[1]]$quality <- quality(reads)[0]
    }

    for (i in seq_along(out)) {
        alignments <- out[[i]]
        read0 <- pattern(alignments)
        adaptor0 <- subject(alignments)
        extended <- .Call(cxx_get_aligned_sequence, aligned(adaptor0), as.character(unaligned(adaptor0)), start(adaptor0), end(adaptor0), aligned(read0))

        output <- DataFrame(score=score(alignments), adaptor=extended[[1]], read=extended[[2]], start=start(read0), end=end(read0))
        if (has.quality) {
            output$quality <- subseq(quality(by.core[[i]]), start=output$start, end=output$end)
        }
        collected[[i]] <- output
    }

    do.call(rbind, collected)
}
