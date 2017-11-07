UMIAlignfastq <- function(UMI, fastqseq, align_pos = 50, score_filter = 8)

#Pairwise alignment UMI-reads
#We start with a pairwise alignment of the reads and the adaptor sequence to get scores and locations of alignments.
#Quality assessment is built in the UMIAlignfastq function and filters by selecting only reads that have an alignment score >8 and in addition an alignment start in the first (fwd) or last (rev) 49 bases. 
#UMIAlignfastq also filters the respective phredscores and stores them in a BStringSet.

{
    UMI_rev <- reverseComplement(UMI)
    align_fwd <- pairwiseAlignment(head(sread(fastqseq)), UMI, type="local-global")
    align_rev <- pairwiseAlignment(head(sread(fastqseq)), UMI_rev, type="local-global")
    
    ind_fwd_pre <- start(pattern(align_fwd))<align_pos&score(align_fwd)>score_filter
    ind_rev_pre <- nchar(head(sread(fastqseq)))-end(pattern(align_rev))<align_pos&score(align_rev)>score_filter
    ind_pre <- ind_fwd_pre!=ind_rev_pre

    ind_fwd <- start(pattern(align_fwd))[ind_pre]<align_pos&score(align_fwd)[ind_pre]>score_filter
    ind_rev <- nchar(head(sread(fastqseq))[ind_pre])-end(pattern(align_rev))[ind_fwd_pre!=ind_rev_pre]<align_pos&score(align_rev)[ind_pre]>score_filter
    
    align_fwd_filt <- align_fwd[which(ind_pre)[ind_fwd]]
    align_rev_filt <- align_rev[which(ind_pre)[ind_rev]]
    
    qual_filt <- head(quality(fastqseq))[ind_pre]
    
    names_filt <- head(ShortRead::id(fastqseq))[ind_pre]
    
    return(list(align_fwd_filt,align_rev_filt, qual_filt, names_filt, ind_pre, ind_fwd, ind_rev))
}
