adaptorAlign <- function(adaptor1, adaptor2, reads, gapOpening=1, gapExtension=5, match=2, mismatch=-5) 
# This function aligns both adaptors to the read sequence with the specified parameters,
# and returns the alignments that best match the sequence (with reverse complementing if necessary).    
{
    if (is.character(adaptor1)) { 
        adaptor1 <- DNAString(adaptor1)
    }
    if (is.character(adaptor2)) {
        adaptor2 <- DNAString(adaptor2)
    }
    adaptor2 <- DNAString(reverse.string(as.character(adaptor2)))
    adaptor1_revcomp <- reverseComplement(adaptor1)
    adaptor2_revcomp <- reverseComplement(adaptor2)
    rev_reads <- reverseComplement(reads)
    
    reads.start <- vector("list", length(reads))
    reads.end <- vector("list", length(reads))
    reads.start.rev <- vector("list", length(reads))
    reads.end.rev <- vector("list", length(reads))
    for(i in 1:length(reads)){
        reads.start[[i]] <- DNAString(str_sub(as.character(reads[[i]]), start = 1, end = 100))
        reads.end[[i]] <- DNAString(str_sub(as.character(reads[[i]]), start = nchar(as.character(reads[[i]]))-100, end = nchar(as.character(reads[[i]]))))
        reads.start.rev[[i]] <- DNAString(str_sub(as.character(rev_reads[[i]]), start = 1, end = 100))
        reads.end.rev[[i]] <- DNAString(str_sub(as.character(rev_reads[[i]]), start = nchar(as.character(rev_reads[[i]]))-100, end = nchar(as.character(rev_reads[[i]]))))
    }
    
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

    # Coercing them to lists.
    align_start <- as.list(align_start)
    align_end <- as.list(align_end)
    align_revcomp_start <- as.list(align_revcomp_start)
    align_revcomp_end <- as.list(align_revcomp_end)

    # Replacing the alignments.
    align_start[is_reverse] <- align_revcomp_end[is_reverse]
    align_end[is_reverse] <- align_revcomp_start[is_reverse]
    reads[is_reverse] <- rev_reads[is_reverse]    

    return(list(adaptor1=align_start, adaptor2=align_end, reads=reads))
}


