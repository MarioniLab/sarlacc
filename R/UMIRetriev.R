UMIRetriev <- function(align_fwd_filt, align_rev_filt, ind_pre, ind_fwd, ind_rev, fastqseq, UMI_fwd_start = 13, UMI_rev_start = 12, UMI_fwd_end = 18, UMI_rev_end = 17)

#Still need to define start and end sequence of N's (UMI)
#Returns forward and reverse (reverse complemented) UMI sequences
    
    #find better way to extract N's!!!
{
    UMI_seq_fwd <- DNAStringSet(str_sub(pattern(align_fwd_filt),UMI_fwd_start, UMI_fwd_end)) #gives out UMI sequence ("NNNNNN") for fwd reads
    UMI_seq_revcomp <- reverseComplement(DNAStringSet(str_sub(pattern(align_rev_filt),UMI_rev_start, UMI_rev_end))) #gives out UMI sequence ("NNNNNN") for rev reads

    names(UMI_seq_fwd) <- head(ShortRead::id(fastqseq))[which(ind_pre)[ind_fwd]]
    names(UMI_seq_revcomp) <- head(ShortRead::id(fastqseq))[which(ind_pre)[ind_rev]]
    
    UMIs <- DNAStringSet(c(UMI_seq_fwd, UMI_seq_revcomp))
    
    return(UMIs)
}
