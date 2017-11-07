chopRead <- function(align_fwd_filt, align_rev_filt, fastqseq, ind_pre, ind_fwd, ind_rev, score_filter = 8)

#Still need to define the score (currently 8)
#Chop read and quality score and return list of both
#After chopping, chopRead reverse complements the rev reads and reverses the phredscores.
#Additionally, all reads and phredscores are named with their corresponding read name from the initial fastq file.
#It return reads and phredscores as a DNAStringSet and BStringSet, respectively (rev reads/phredscores have been added to the StringSet behind the fwd reads/phredscores)

{
    start_fwd <- start(pattern(align_fwd_filt[score(align_fwd_filt)>score_filter]))
    end_fwd <- end(pattern(align_fwd_filt[score(align_fwd_filt)>score_filter]))
    start_rev <- start(pattern(align_rev_filt[score(align_rev_filt)>score_filter]))
    end_rev <- end(pattern(align_rev_filt[score(align_rev_filt)>score_filter]))
        
    chop_read_fwd <- DNAStringSet(str_sub(head(sread(fastqseq))[which(ind_pre)[ind_fwd]], end_fwd+1))
    chop_read_rev <- DNAStringSet(str_sub(head(sread(fastqseq))[which(ind_pre)[ind_rev]], end=start_rev-1))
        
    chop_qual_fwd <- BStringSet(str_sub(quality(getClass(head(quality(fastqseq))))[which(ind_pre)[ind_fwd]], end_fwd+1))
    chop_qual_rev <- BStringSet(str_sub(quality(getClass(head(quality(fastqseq))))[which(ind_pre)[ind_rev]], end=start_rev-1))
    
    
    #Final sorting of reads with names and construct UMI and read list and reverse the rev reads and the respective quality scores

    #reverse reads and quality scores
    chop_read_rev_rev <- reverseComplement(chop_read_rev)
    chop_qual_rev_rev <- BStringSet(reverse.string(as.character(chop_qual_rev)))

    names(chop_read_fwd) <- head(ShortRead::id(fastqseq))[which(ind_pre)[ind_fwd]]
    names(chop_read_rev_rev) <- head(ShortRead::id(fastqseq))[which(ind_pre)[ind_rev]]
    names(chop_qual_fwd) <- head(ShortRead::id(fastqseq))[which(ind_pre)[ind_fwd]]
    names(chop_qual_rev_rev) <- head(ShortRead::id(fastqseq))[which(ind_pre)[ind_rev]]

    reads <- DNAStringSet(c(chop_read_fwd, chop_read_rev_rev))
    qualities <- c(chop_qual_fwd, chop_qual_rev_rev)
    
    writeXStringSet(reads, "/Users/florian/Desktop/Cambridge_Bioinformatics/Pipeline_project/chopped_readsq.fastq", format = "fastq", qualities= qualities) 

    return(list(reads, qualities))
}
