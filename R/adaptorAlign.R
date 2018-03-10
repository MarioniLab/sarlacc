#' @export
adaptorAlign <- function(adaptor1, adaptor2, reads, tolerance=100, gapOpening=1, gapExtension=5, match=1, mismatch=0)
# This function aligns both adaptors to the read sequence with the specified parameters,
# and returns the alignments that best match the sequence (with reverse complementing if necessary).
#    
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 7 November 2017    
{
    pre.out <- .preprocess_input(adaptor1, adaptor2, reads)
    adaptor1 <- pre.out$adaptor1
    adaptor2 <- pre.out$adaptor2
    reads <- pre.out$reads

    reads.out <- .get_front_and_back(reads, tolerance)
    reads.start <- reads.out$front
    reads.end <- reads.out$back

    .internal_align(adaptor1, adaptor2, reads.start, reads.end,
                    gapOpening=gapOpening, gapExtension=gapExtension,
                    match=match, mismatch=mismatch)
}

#' @importFrom Biostrings DNAString 
.preprocess_input <- function(adaptor1, adaptor2, reads) 
# Coerces all inputs to DNAString or DNAStringSet objects.
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
    if (is.null(names(reads))) { 
        names(reads) <- paste0("READ", seq_along(reads))
    }
    return(list(adaptor1=adaptor1, adaptor2=adaptor2, reads=reads))
}

#' @importFrom Biostrings reverseComplement
#' @importFrom XVector subseq
.get_front_and_back <- function(reads, tolerance) 
# Extracting the front part of the read (forward strand) 
# and the end of the read (reverse strand).
{
    tolerance <- pmin(tolerance, width(reads))
    reads.start <- subseq(reads, start = 1, width = tolerance)
    reads.end <- subseq(reads, end = width(reads), width = tolerance)
    reads.end <- reverseComplement(reads.end)
    return(list(front=reads.start, back=reads.end))
}

#' @importFrom Biostrings DNAStringSet reverseComplement 
#' @importFrom BiocGenerics score
#' @importFrom methods is
.internal_align <- function(adaptor1, adaptor2, reads.start, reads.end, gapOpening=1, gapExtension=5, match=1, mismatch=0) 
# Takes adaptor sequences on the forward strand, the start of the read
# on the forward strand, and the end of the read on the reverse strand,
# and performs alignments of the adaptors to the read start/ends.
{
    has.quality <- is(reads, "QualityScaledDNAStringSet")
    all.args <- .setup_alignment_args(has.quality, gapOpening, gapExtension, match, mismatch)
    align.out <- .get_all_alignments(adaptor1, adaptor2, reads.start, reads.end, all.args=all.args)
    align_start <- align.out$start
    align_end <- align.out$end
    align_revcomp_start <- align.out$rc.start
    align_revcomp_end <- align.out$rc.end
    
    # Figuring out the strand.
    fscore <- pmax(score(align_start), score(align_end))
    rscore <- pmax(score(align_revcomp_start), score(align_revcomp_end))
    is_reverse <- fscore < rscore

    # Extracting all the useful information out of them.
    if (has.quality) {
        qual.start <- quality(reads.start)
        qual.end <- quality(reads.end)
    } else {
        qual.start <- NULL
        qual.end <- NULL
    }

    align_start <- .align_info_extractor(align_start, quality=qual.start)
    align_end <- .align_info_extractor(align_end, quality=qual.end)
    align_revcomp_start <- .align_info_extractor(align_revcomp_start, quality=qual.end)
    align_revcomp_end <- .align_info_extractor(align_revcomp_end, quality=qual.start)

    # Replacing the alignments.
    align_start[is_reverse,] <- align_revcomp_start[is_reverse,]
    align_end[is_reverse,] <- align_revcomp_end[is_reverse,]
    reads[is_reverse] <- reverseComplement(reads[is_reverse])

    # Adjusting the reverse coordinates for the read length.
    old.start <- align_end$start
    old.end <- align_end$end
    align_end$start <- width(reads) - old.start + 1L
    align_end$end <- width(reads) - old.end + 1L
       
    rownames(align_start) <- rownames(align_end) <- names(reads) 
    names(is_reverse) <- names(reads)
    return(list(adaptor1=align_start, adaptor2=align_end, reads=reads, reversed=is_reverse))
}

#' @importFrom Biostrings nucleotideSubstitutionMatrix
.setup_alignment_args <- function(has.quality, gapOpening, gapExtension, match, mismatch) {
    all.args <- list(type="local-global", gapOpening=gapOpening, gapExtension=gapExtension)
    if (!has.quality) { 
        all.args$substitutionMatrix <- nucleotideSubstitutionMatrix(match=match, mismatch=mismatch)
    }
    return(all.args)
}

#' @importFrom Biostrings pairwiseAlignment
.get_all_alignments <- function(adaptor1, adaptor2, reads.start, reads.end, all.args=list(), ...) {
    align_start <- do.call(pairwiseAlignment, c(list(pattern=reads.start, subject=adaptor1, ...), all.args))
    align_end <- do.call(pairwiseAlignment, c(list(pattern=reads.end, subject=adaptor2, ...), all.args))
    align_revcomp_start <- do.call(pairwiseAlignment, c(list(pattern=reads.end, subject=adaptor1, ...), all.args))
    align_revcomp_end <- do.call(pairwiseAlignment, c(list(pattern=reads.start, subject=adaptor2, ...), all.args))
    list(start=align_start, end=align_end, rc.start=align_revcomp_start, rc.end=align_revcomp_end)
}

#' @importFrom Biostrings pattern subject
#' @importFrom BiocGenerics score start end
#' @importFrom S4Vectors DataFrame
#' @importFrom XVector subseq
.align_info_extractor <- function(alignments, quality=NULL) {
    P <- pattern(alignments)
    S <- subject(alignments)
    output <- DataFrame(score=score(alignments), adaptor=as.character(P), 
                        read=as.character(S), start=start(P), end=end(P))
    if (!is.null(quality)) {
        pattern.qual <- subseq(quality, start=output$start, end=output$end)
        output$quality <- pattern.qual
    }    
    return(output)
}
