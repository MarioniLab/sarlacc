adaptorAlign <- function(adaptor1, adaptor2, reads, quality = NULL, gapOpening=1, gapExtension=5, match=2, mismatch=-5, tolerance=100) 
# This function aligns both adaptors to the read sequence with the specified parameters,
# and returns the alignments that best match the sequence (with reverse complementing if necessary).
{
    # Checking the sensibility of the adaptors.
    if (is.character(adaptor1)) { 
        adaptor1 <- DNAString(adaptor1)
    }
    if (is.character(adaptor2)) {
        adaptor2 <- DNAString(adaptor2)
    }

    # Checking the sensibility of the reads.
    if (is.character(reads)) {
        reads <- DNAString(reads)
    }
    if (is.null(names(reads))) { 
        names(reads) <- paste0("READ", seq_along(reads))
    }     
   
    tolerance <- pmin(tolerance, width(reads))
    reads.start <- subseq(reads, start = 1, width = tolerance)
    reads.end <- subseq(reads, end = width(reads), width = tolerance)
    reads.end <- reverseComplement(reads.end)
    
    # Aligning all sequences.
    submat <- nucleotideSubstitutionMatrix(match=match, mismatch=mismatch)
    all.args <- list(type="local-global", gapOpening=gapOpening, gapExtension=gapExtension, substitutionMatrix=submat)
    align_start <- do.call(pairwiseAlignment, c(list(pattern=reads.start, subject=adaptor1), all.args))
    align_end <- do.call(pairwiseAlignment, c(list(pattern=reads.end, subject=adaptor2), all.args))
    align_revcomp_start <- do.call(pairwiseAlignment, c(list(pattern=reads.end, subject=adaptor1), all.args))
    align_revcomp_end <- do.call(pairwiseAlignment, c(list(pattern=reads.start, subject=adaptor2), all.args))
    
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
    align_start[is_reverse,] <- align_revcomp_start[is_reverse,]
    align_end[is_reverse,] <- align_revcomp_end[is_reverse,]
    reads[is_reverse] <- reverseComplement(reads[is_reverse])

    # Adjusting the reverse coordinates for the read length.
    old.start <- align_end$start.pattern
    old.end <- align_end$end.pattern
    align_end$start.pattern <- width(reads) - old.start + 1L
    align_end$end.pattern <- width(reads) - old.end + 1L
       
    # Read quality, if provided.
    if (!is.null(quality)){
        if (is.character(quality)) { 
            quality <- BStringSet(quality)
        }
        quality[is_reverse] <- reverse(quality[is_reverse])
        names(quality) <- names(reads)
    }
   
    rownames(align_start) <- rownames(align_end) <- names(reads) 
    names(is_reverse) <- names(reads)
    return(list(adaptor1=align_start, adaptor2=align_end, reads=reads, quality=quality, reversed=is_reverse))
}

.align_info_extractor <- function(alignments) {
    P <- pattern(alignments)
    S <- subject(alignments)
    DataFrame(pattern=as.character(P), start.pattern=start(P), end.pattern=end(P),
              subject=as.character(S), start.subject=start(S), end.subject=end(S),
              score=score(alignments))
}
