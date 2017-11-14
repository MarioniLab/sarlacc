UMIRetriev <- function(align_start_filt, align_end_filt = NULL, align_revcomp_start_filt = NULL, align_revcomp_end_filt, ind, ind_revcomp, fastqseq, UMI_pos1 = c(20, 31), UMI_pos2 = c(12, 24), UMI_pos3 = NULL, UMI_pos4 = NULL)
    #Still need to define start and end sequence of N's (UMI)
    #Returns forward and reverse (reverse complemented) UMI sequences bound behind each other (fwd - rev) in DNAStringSet
{
    
    if((is.null(UMI_pos3)) | (as.character(adapt1) == as.character(reverseComplement(adapt2)))){
        
        UMI_seq_fwd <- .UMIindelextract(alignment = align_start_filt, UMI_pos = UMI_pos1)
        UMI_seq_rev <- .UMIindelextract(alignment = align_revcomp_end_filt, UMI_pos = UMI_pos2)
        
        UMIfwd.collide <- .UMIcollide(UMIindelextract = UMI_seq_fwd)
        UMIrev.collide <- .UMIcollide(UMIindelextract = UMI_seq_rev)
        
        names(UMIfwd.collide) <- ShortRead::id(fastqseq)[ind]
        names(UMIrev.collide) <- ShortRead::id(fastqseq)[ind_revcomp]
        
        UMIs <- DNAStringSet(as.character(c(UMIfwd.collide, UMIrev.collide)))
        
    }else{
        
        UMI_seq_fwd <- .UMIindelextract(alignment = align_start_filt, UMI_pos = UMI_pos1)
        UMI_seq_rev <- .UMIindelextract(alignment = align_revcomp_end_filt, UMI_pos = UMI_pos2)
        UMI_seq_2_fwd <- .UMIindelextract(alignment = align_end_filt, UMI_pos = UMI_pos3)
        UMI_seq_2_rev <- .UMIindelextract(alignment = align_revcomp_start_filt, UMI_pos = UMI_pos4)
        
        UMIf.collide <- .UMIcollide(UMIindelextract = UMI_seq_fwd)
        UMIr.collide <- .UMIcollide(UMIindelextract = UMI_seq_rev)
        UMI2f.collide <- .UMIcollide(UMIindelextract = UMI_seq_2_fwd)
        UMI2r.collide <- .UMIcollide(UMIindelextract = UMI_seq_2_rev)
        
        names(UMIf.collide) <- ShortRead::id(fastqseq)[ind]
        names(UMIr.collide) <- ShortRead::id(fastqseq)[ind_revcomp]
        names(UMI2f.collide) <- ShortRead::id(fastqseq)[ind]
        names(UMI2r.collide) <- ShortRead::id(fastqseq)[ind_revcomp]
        
        UMIs <- DNAStringSet(as.character(c(UMIf.collide, UMIr.collide, UMI2f.collide, UMI2r.collide)))
    }
    
    .UMIindelextract <- function(alignment, UMI_pos)
    {
        UMI_seq <- vector("list", length(alignment))
        for (i in 1:length(alignment)){
            UMI_fwd_start1 <- UMI_pos[1] + sum(width(indel(subject(alignment)))[i][start(indel(subject(alignment)))[i] <= UMI_pos[1]])
            UMI_fwd_end1 <- UMI_pos[2] + sum(width(indel(subject(alignment)))[i][start(indel(subject(alignment)))[i] <= UMI_pos[1]])
            UMI_seq[[i]] <- DNAStringSet(str_sub(pattern(alignment[i]),UMI_fwd_start1, UMI_fwd_end1))
        }
        
        return(UMI_seq)
    }
    
    .UMIcollide <- function(UMIindelextract)
    {
        UMI.collide <- vector("list", length(UMIindelextract))
        for (i in 1:length(UMI.collide)){
            cuts <- grep("[A-Z]", strsplit(as.character(UMIindelextract[[i]]), NULL)[[1]])
            x <- vector("list", length(cuts))
            for (a in 1:length(cuts)){
                x[a] <- strsplit(as.character(UMIindelextract[[i]]), NULL)[[1]][cuts[a]]
            }
            UMI.collide[i] <- paste(x, collapse = "")
        } 
        
        return(UMI.collide)
    }
    
    
    return(UMIs)
}

