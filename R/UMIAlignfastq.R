UMIAlignfastq <- function(adaptor1, adaptor2, fastqseq, align_pos = 40, score_filter = -12)
    #For short adaptor1's (15, 6) set score_filter to 8 
    #For long adaptor1's (45, 12) set score_filter to -12
    #adaptors must be entered in 5'->3' direction

    {
        
        adaptor1_revcomp <- reverseComplement(adaptor1)
        adaptor2_revcomp <- reverseComplement(adaptor2)
        
        align_start <- pairwiseAlignment(sread(fastqseq), adaptor1, type="local-global")
        align_end <- pairwiseAlignment(sread(fastqseq), adaptor2, type="local-global")
        align_revcomp_end <- pairwiseAlignment(sread(fastqseq), adaptor1_revcomp, type="local-global")
        align_revcomp_start <- pairwiseAlignment(sread(fastqseq), adaptor2_revcomp, type="local-global")
        
        ind <- (start(pattern(align_start)) < align_pos & score(align_start) > score_filter) & (nchar(sread(fastqseq))-end(pattern(align_end)) < align_pos & score(align_end) > score_filter)
        ind_revcomp <- (start(pattern(align_revcomp_start)) < align_pos & score(align_revcomp_start) > score_filter) & (nchar(sread(fastqseq))-end(pattern(align_revcomp_end)) < align_pos & score(align_revcomp_end) > score_filter)
        
        if(sum(ind==ind_revcomp)==length(sread(fastqseq))){
            ind_comb <- ind
        }else{
            ind_comb <- ind!=ind_revcomp
            ind[!ind_comb] <- FALSE
            ind_revcomb[!ind_comb] <- FALSE
        }
        
        align_start_filt <- align_start[ind]
        align_end_filt <- align_end[ind]
        align_revcomp_start_filt <- align_revcomp_start[ind_revcomp] 
        align_revcomp_end_filt <- align_revcomp_end[ind_revcomp] 
        
        
        qual_filt <- quality(fastqseq)[ind_comb]
        
        names_filt <- ShortRead::id(fastqseq)[ind_comb]
        
        return(list(align_start_filt, align_end_filt, align_revcomp_start_filt, align_revcomp_end_filt, qual_filt, names_filt, ind_comb, ind, ind_revcomp))
}
