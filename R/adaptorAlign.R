adaptorAlign <- function(adaptor1, adaptor2, reads, quality = NULL, gapOpening=1, gapExtension=5, match=2, mismatch=-5, tolerance=100) 
# This function aligns both adaptors to the read sequence with the specified parameters,
# and returns the alignments that best match the sequence (with reverse complementing if necessary).    
{
    if (is.character(adaptor1)) { 
        adaptor1 <- DNAString(adaptor1)
    }
    if (is.character(adaptor2)) {
        adaptor2 <- DNAString(adaptor2)
    }
    adaptor1_revcomp <- reverseComplement(adaptor1)
    adaptor2_revcomp <- reverseComplement(adaptor2)

    tolerance <- pmin(tolerance, width(reads))
    reads.start <- subseq(reads, start = 1, width = tolerance)
    reads.end <- subseq(reads, end = width(reads), width = tolerance)
    reads.start.rev <- reverseComplement(reads.end)
    reads.end.rev <- reverseComplement(reads.start)
    
    # Aligning all sequences.
    submat <- nucleotideSubstitutionMatrix(match=match, mismatch=mismatch)
    all.args <- list(type="local-global", gapOpening=gapOpening, gapExtension=gapExtension, substitutionMatrix=submat)
    align_start <- do.call(pairwiseAlignment, c(list(pattern=reads.start, subject=adaptor1), all.args))
    align_end <- do.call(pairwiseAlignment, c(list(pattern=reads.end, subject=adaptor2_revcomp), all.args))
    align_revcomp_start <- do.call(pairwiseAlignment, c(list(pattern=reads.start.rev, subject=adaptor2), all.args))
    align_revcomp_end <- do.call(pairwiseAlignment, c(list(pattern=reads.end.rev, subject=adaptor1_revcomp), all.args))
    
    # Figuring out the strand.
    fscore <- pmax(score(align_start), score(align_end))
    rscore <- pmax(score(align_revcomp_start), score(align_revcomp_end))
    is_reverse <- fscore < rscore

    # Extracting all the useful information out of them.
    align_start <- .align_info_extractor(align_start)
    align_end <- .align_info_extractor(align_end)
    align_revcomp_start <- .align_info_extractor(align_revcomp_start)
    align_revcomp_end <- .align_info_extractor(align_revcomp_end)

    # Replacing the alignments.
    align_start[is_reverse,] <- align_revcomp_end[is_reverse,]
    align_end[is_reverse,] <- align_revcomp_start[is_reverse,]
    reads[is_reverse] <- reverseComplement(reads[is_reverse])

    # Adjusting the reverse coordinates for the read length.
    adjustments <- width(reads) - tolerance # guaranteed to be at least zero, by pmin() above.
    align_end$start.pattern <- align_end$start.pattern + adjustments
    align_end$end.pattern <- align_end$end.pattern + adjustments
       
    # Read quality, if provided.
    if (!is.null(quality)){
        if (is.character(quality)) { 
            quality <- BStringSet(quality)
        }
        quality[is_reverse] <- reverse(quality[is_reverse])
    }
    
    return(list(adaptor1=align_start, adaptor2=align_end, reads=reads, quality=quality, reversed=is_reverse))
}

.align_info_extractor <- function(alignments) {
    P <- pattern(alignments)
    S <- subject(alignments)
    DataFrame(pattern=as.character(P), start.pattern=start(P), end.pattern=end(P),
              subject=as.character(S), start.subject=start(S), end.subject=end(S),
              score=score(alignments))
}
