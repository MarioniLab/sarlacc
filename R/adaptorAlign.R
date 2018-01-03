adaptorAlign <- function(adaptor1, adaptor2, reads, gapOpening=1, gapExtension=5, match=2, mismatch=-5, tolerance=100) 
# This function aligns both adaptors to the read sequence with the specified parameters,
# and returns the alignments that best match the sequence (with reverse complementing if necessary).
#    
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 7 November 2017    
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

    # Taking subsequences of the reads for pairwise alignment.
    tolerance <- pmin(tolerance, width(reads))
    reads.start <- subseq(reads, start = 1, width = tolerance)
    reads.end <- subseq(reads, end = width(reads), width = tolerance)
    reads.end <- reverseComplement(reads.end)
    
    # Aligning all sequences to the adaptor.
    all.args <- list(type="local-global", gapOpening=gapOpening, gapExtension=gapExtension)
    if (!is(reads, "QualityScaledDNAStringSet")) { 
        all.args$substitutionMatrix <- nucleotideSubstitutionMatrix(match=match, mismatch=mismatch)
    }

    align_start <- do.call(pairwiseAlignment, c(list(pattern=reads.start, subject=adaptor1), all.args))
    align_end <- do.call(pairwiseAlignment, c(list(pattern=reads.end, subject=adaptor2), all.args))
    align_revcomp_start <- do.call(pairwiseAlignment, c(list(pattern=reads.end, subject=adaptor1), all.args))
    align_revcomp_end <- do.call(pairwiseAlignment, c(list(pattern=reads.start, subject=adaptor2), all.args))
    
    # Figuring out the strand.
    fscore <- pmax(score(align_start), score(align_end))
    rscore <- pmax(score(align_revcomp_start), score(align_revcomp_end))
    is_reverse <- fscore < rscore

    # Extracting all the useful information out of them.
    align_start <- .align_info_extractor(align_start, quality=quality(reads.start))
    align_end <- .align_info_extractor(align_end, quality=quality(reads.end))
    align_revcomp_start <- .align_info_extractor(align_revcomp_start, quality=quality(reads.end))
    align_revcomp_end <- .align_info_extractor(align_revcomp_end, quality=quality(reads.start))

    # Replacing the alignments.
    align_start[is_reverse,] <- align_revcomp_start[is_reverse,]
    align_end[is_reverse,] <- align_revcomp_end[is_reverse,]
    reads[is_reverse] <- reverseComplement(reads[is_reverse])

    # Adjusting the reverse coordinates for the read length.
    old.start <- align_end$start.pattern
    old.end <- align_end$end.pattern
    align_end$start.pattern <- width(reads) - old.start + 1L
    align_end$end.pattern <- width(reads) - old.end + 1L
       
    rownames(align_start) <- rownames(align_end) <- names(reads) 
    names(is_reverse) <- names(reads)
    return(list(adaptor1=align_start, adaptor2=align_end, reads=reads, reversed=is_reverse))
}

.align_info_extractor <- function(alignments, quality=NULL) {
    P <- pattern(alignments)
    S <- subject(alignments)
    output <- DataFrame(score=score(alignments), adaptor=as.character(P), 
                        read=as.character(S), start=start(P), end=end(P))
    if (!is.null(quality)) {
        pattern.qual <- subseq(quality, start=output$start.pattern, end=output$end.pattern)
        output$quality <- pattern.qual
    }    
    return(output)
}
