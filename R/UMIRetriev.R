UMIRetriev <- function(align_fwd_filt, align_rev_filt, ind_fwd, ind_rev, fastqseq, UMI_fwd_start = 13, UMI_rev_start = 12, UMI_fwd_end = 18, UMI_rev_end = 17)
    #Still need to define start and end sequence of N's (UMI)
    #Returns forward and reverse (reverse complemented) UMI sequences bound behind each other (fwd - rev) in DNAStringSet
{
    UMI_seq_fwd <- vector("list", length(align_fwd_filt))
    for (i in 1:length(align_fwd_filt)){
        if (start(indel(subject(align_fwd_filt)))[[i]] <= UMI_fwd_start){
            UMI_fwd_start1 <- UMI_fwd_start + width(indel(subject(align_fwd_filt)))
            UMI_fwd_end1 <- UMI_fwd_end + width(indel(subject(align_fwd_filt)))
            UMI_seq_fwd[i] <- DNAStringSet(str_sub(pattern(align_fwd_filt),UMI_fwd_start1, UMI_fwd_end1))
        } else{
            UMI_seq_fwd[i] <- DNAStringSet(str_sub(pattern(align_fwd_filt),UMI_fwd_start, UMI_fwd_end))
        }
    }
    
    UMI_seq_revcomp <- vector("list", length(align_rev_filt))
    for (i in 1:length(align_rev_filt)){
        if (start(indel(subject(align_rev_filt)))[[i]] >= UMI_rev_start){
            UMI_rev_start1 <- UMI_rev_start + width(indel(subject(align_rev_filt)))
            UMI_rev_end1 <- UMI_rev_end + width(indel(subject(align_rev_filt)))
            UMI_seq_revcomp[i] <- reverseComplement(DNAStringSet(str_sub(pattern(align_rev_filt),UMI_rev_start1, UMI_rev_end1)))
        } else{
            UMI_seq_revcomp[i] <- reverseComplement(DNAStringSet(str_sub(pattern(align_rev_filt),UMI_rev_start, UMI_rev_end)))
        }
    }
    
    UMIf.collide <- vector("list", length(UMI_seq_fwd))
    for (i in 1:length(UMI_seq_fwd)){
        cuts <- grep("[A-Z]", strsplit(as.character(UMI_seq_fwd[[i]]), NULL)[[i]])
        x <- vector("list", )
        for (a in 1:length(cuts)){
            x[a] <- strsplit(as.character(UMI_seq_fwd[[i]]), NULL)[[i]][cuts[a]]
        }
        UMIf.collide[i] <- paste(x, collapse = "")
    }
    
    UMIr.collide <- vector("list", length(UMI_seq_revcomp))
    for (i in 1:length(UMI_rev.collide)){
        cuts <- grep("[A-Z]", strsplit(as.character(UMI_seq_revcomp[[i]]), NULL)[[i]])
        x <- vector("list", )
        for (a in 1:length(cuts)){
            x[a] <- strsplit(as.character(UMI_seq_revcomp[[i]]), NULL)[[i]][cuts[a]]
        }
        UMIr.collide[i] <- paste(x, collapse = "")
    }
    
    
    names(UMI_fwd.collide) <- head(ShortRead::id(fastqseq))[ind_fwd]
    names(UMI_rev.collide) <- head(ShortRead::id(fastqseq))[ind_rev]
    
    UMIs <- DNAStringSet(as.character(c(UMIf.collide, UMIr.collide)))
    
    
    return(UMIs)
}