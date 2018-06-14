#' @export
#' @importFrom Biostrings reverseComplement 
#' @importFrom BiocGenerics score
#' @importFrom S4Vectors metadata
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom methods is
adaptorAlign <- function(adaptor1, adaptor2, reads, tolerance=100, gapOpening=1, gapExtension=5, match=1, mismatch=0)
# This function aligns both adaptors to the read sequence with the specified parameters,
# and returns the alignments that best match the sequence (with reverse complementing if necessary).
#    
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 7 November 2017    
{
    pre.out <- .preprocess_input(adaptor1, adaptor2, reads, add.names=TRUE)
    adaptor1 <- pre.out$adaptor1
    adaptor2 <- pre.out$adaptor2
    reads <- pre.out$reads

    # Getting the start and (rc'd) end of the read.
    reads.out <- .get_front_and_back(reads, tolerance)
    reads.start <- reads.out$front
    reads.end <- reads.out$back

    # Performing the alignment of each adaptor to the start/end of the read.
    has.quality <- is(reads, "QualityScaledDNAStringSet")
    all.args <- .setup_alignment_args(has.quality, gapOpening, gapExtension, match, mismatch)
    align.out <- .get_all_alignments(adaptor1, adaptor2, reads.start, reads.end, all.args=all.args)
    align_start <- align.out$start
    align_end <- align.out$end
    align_revcomp_start <- align.out$rc.start
    align_revcomp_end <- align.out$rc.end
    
    # Figuring out the strand.
    strand.out <- .resolve_strand(score(align_start), score(align_end), score(align_revcomp_start), score(align_revcomp_end))
    is_reverse <- strand.out$reversed

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

    metadata(align_start)$sequence <- adaptor1
    metadata(align_end)$sequence <- adaptor2

    # Adjusting the reverse coordinates for the read length.
    old.start <- align_end$start
    old.end <- align_end$end
    align_end$start <- width(reads) - old.start + 1L
    align_end$end <- width(reads) - old.end + 1L
       
    rownames(align_start) <- rownames(align_end) <- names(reads) 
    names(is_reverse) <- names(reads)
    return(list(adaptor1=align_start, adaptor2=align_end, reads=reads, reversed=is_reverse,
                parameters=list(tolerance=tolerance, gapOpening=gapOpening, gapExtension=gapExtension, 
                                match=match, mismatch=mismatch)))
}

#' @importFrom Biostrings DNAString DNAStringSet
.preprocess_input <- function(adaptor1, adaptor2, reads, add.names=FALSE) 
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
    if (add.names && is.null(names(reads))) { 
        names(reads) <- paste0("READ", seq_along(reads))
    }
    return(list(adaptor1=adaptor1, adaptor2=adaptor2, reads=reads))
}

#' @importFrom Biostrings reverseComplement
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

#' @importFrom Biostrings nucleotideSubstitutionMatrix
.setup_alignment_args <- function(has.quality, gapOpening, gapExtension, match, mismatch) {
    all.args <- list(type="local-global", gapOpening=gapOpening, gapExtension=gapExtension)
    if (!has.quality) { 
        all.args$substitutionMatrix <- nucleotideSubstitutionMatrix(match=match, mismatch=mismatch)
    }
    return(all.args)
}

#' @importFrom Biostrings pairwiseAlignment
.get_all_alignments <- function(adaptor1, adaptor2, reads.start, reads.end, all.args=list(), ...) 
# This performs the alignment of adaptor 1 to the start (start) and adaptor 2 to the end (end),
# or adaptor 1 to the end (rc.start, as it's the start of the reverse complemented read) and 
# adaptor 2 to the start (rc.end, as it's the end of the reverse complemented read).
{
    align_start <- do.call(pairwiseAlignment, c(list(pattern=reads.start, subject=adaptor1, ...), all.args))
    align_end <- do.call(pairwiseAlignment, c(list(pattern=reads.end, subject=adaptor2, ...), all.args))
    align_revcomp_start <- do.call(pairwiseAlignment, c(list(pattern=reads.end, subject=adaptor1, ...), all.args))
    align_revcomp_end <- do.call(pairwiseAlignment, c(list(pattern=reads.start, subject=adaptor2, ...), all.args))
    list(start=align_start, end=align_end, rc.start=align_revcomp_start, rc.end=align_revcomp_end)
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

#' @importFrom Biostrings pattern subject aligned unaligned
#' @importFrom BiocGenerics score 
#' @importFrom stats start end
#' @importFrom S4Vectors DataFrame
#' @importFrom XVector subseq
.align_info_extractor <- function(alignments, quality=NULL) {
    read0 <- pattern(alignments)
    adaptor0 <- subject(alignments)
    extended <- .Call(cxx_get_aligned_sequence, aligned(adaptor0), as.character(unaligned(adaptor0)), start(adaptor0), end(adaptor0), aligned(read0))

    output <- DataFrame(score=score(alignments), adaptor=extended[[1]], read=extended[[2]], start=start(read0), end=end(read0))
    if (!is.null(quality)) {
        pattern.qual <- subseq(quality, start=output$start, end=output$end)
        output$quality <- pattern.qual
    }    
    return(output)
}


