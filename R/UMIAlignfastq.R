UMIAlignfastq <- function(UMI, fastqseq, align_pos = 50, score_filter = 8)
    
    #Pairwise alignment UMI-reads
    #We start with a pairwise alignment of the reads and the adaptor sequence to get scores and locations of alignments.
    #Quality assessment is built in the UMIAlignfastq function and filters by selecting only reads that have an alignment score >8 and in addition an alignment start in the first (fwd) or last (rev) 49 bases. 
    #UMIAlignfastq also filters the respective phredscores and stores them in a BStringSet.
    
{
    UMI_rev <- reverseComplement(UMI)
    align_fwd <- pairwiseAlignment(sread(fastqseq), UMI, type="local-global")
    align_rev <- pairwiseAlignment(sread(fastqseq), UMI_rev, type="local-global")
    
    ind_fwd <- start(pattern(align_fwd)) < align_pos & score(align_fwd) > score_filter
    ind_rev <- nchar(sread(fastqseq))-end(pattern(align_rev)) < align_pos&score(align_rev) > score_filter
    ind_pre <- ind_fwd!=ind_rev
    
    ind_fwd[!ind_pre] <- FALSE
    ind_rev[!ind_pre] <- FALSE
    
    align_fwd_filt <- align_fwd[ind_fwd]
    align_rev_filt <- align_rev[ind_rev]
    
    qual_filt <- quality(fastqseq)[ind_pre]
    
    names_filt <- ShortRead::id(fastqseq)[ind_pre]
    
    return(list(align_fwd_filt,align_rev_filt, qual_filt, names_filt, ind_pre, ind_fwd, ind_rev))
}
