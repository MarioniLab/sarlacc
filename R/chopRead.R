chopRead <- function(align_start_filt, align_end_filt, align_revcomp_start_filt, align_revcomp_end_filt, fastqseq, ind_comb, ind, ind_revcomp, score_filter = 8)
    
    #Still need to define the score (currently 8)
    #Chop read and quality score and return list of both
    #After chopping, chopRead reverse complements the rev reads and reverses the phredscores.
    #Additionally, all reads and phredscores are named with their corresponding read name from the initial fastq file.
    #It return reads and phredscores as a DNAStringSet and BStringSet, respectively (rev reads/phredscores have been added to the StringSet behind the fwd reads/phredscores)
    
{
    starta1 <- start(pattern(align_start_filt[score(align_start_filt)>score_filter]))
    enda1 <- end(pattern(align_start_filt[score(align_start_filt)>score_filter]))
    starta2 <- start(pattern(align_end_filt[score(align_end_filt)>score_filter]))
    enda2 <- end(pattern(align_end_filt[score(align_end_filt)>score_filter]))
    
    starta1_revcomp <- start(pattern(align_revcomp_start_filt[score(align_revcomp_start_filt)>score_filter]))
    enda1_revcomp <- end(pattern(align_revcomp_start_filt[score(align_revcomp_start_filt)>score_filter]))
    starta2_revcomp<- start(pattern(align_revcomp_end_filt[score(align_revcomp_end_filt)>score_filter]))
    enda2_revcomp <- end(pattern(align_revcomp_end_filt[score(align_revcomp_end_filt)>score_filter]))
    
    
    chop_read <- DNAStringSet(str_sub(fastqseq[ind], start = enda1+1, end = starta2-1))
    chop_read_revcomp <- DNAStringSet(str_sub(fastqseq[ind_revcomp], start = enda1_revcomp+1, end = starta2_revcomp-1))
    
    chop_qual <- BStringSet(str_sub(quality(getClass(quality(sread(fastqseq))))[ind], start = enda1+1, end = starta2-1))
    chop_qual_revcomp <- BStringSet(str_sub(quality(getClass(quality(sread(fastqseq))))[ind_revcomp], start = enda1_revcomp+1, end = starta2_revcomp-1))
    
    
    #Final sorting of reads with names and construct UMI and read list and reverse the rev reads and the respective quality scores
    
    #reverse reads and quality scores
    
    if(as.character(adaptor1) == as.character(reverseComplement(adaptor2))){
        names(chop_read) <- ShortRead::id(fastqseq)[ind]
        names(chop_read_revcomp) <- ShortRead::id(fastqseq)[ind_revcomp]
        names(chop_qual) <- ShortRead::id(fastqseq)[ind]
        names(chop_qual_revcomp) <- ShortRead::id(fastqseq)[ind_revcomp]
        
        reads <- chop_read
        qualities <- chop_qual
    }else{
        chop_read_rev_rev <- reverseComplement(chop_read_revcomp)
        chop_qual_rev_rev <- BStringSet(reverse.string(as.character(chop_qual_revcomp)))
        
        names(chop_read) <- ShortRead::id(fastqseq)[ind]
        names(chop_read_rev_rev) <- ShortRead::id(fastqseq)[ind_revcomp]
        names(chop_qual) <- ShortRead::id(fastqseq)[ind]
        names(chop_qual_rev_rev) <- ShortRead::id(fastqseq)[ind_revcomp]
        
        reads <- DNAStringSet(c(chop_read, chop_read_rev_rev))
        qualities <- c(chop_qual, chop_qual_rev_rev)
    }
    
    return(list(reads, qualities))
}