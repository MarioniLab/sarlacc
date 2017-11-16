UMIRetriev <- function(chop_data, UMI_adaptor1_pos = NULL, UMI_adaptor2_pos = NULL)

{
    if(is.null(UMI_adaptor1_pos)){
        
        UMI1.collide <- NULL
        
    }else{
        
        UMI_seq1 <- .UMIindelextract(alignment = chop_data$adaptor1filt, UMI_pos = UMI_adaptor1_pos)
        UMI1.collide <- .UMIcollide(UMIseq = UMI_seq1)
        
        names(UMI1.collide) <- names(chop_data$chopread)
        
        UMI1 <- DNAStringSet(as.character(UMI1.collide))
    }
    
    if(is.null(UMI_adaptor2_pos)){
        
        UMI2.collide <- NULL
        
    }else{
        UMI_seq2 <- .UMIindelextract(alignment = chop_data$adaptor2filt, UMI_pos = UMI_adaptor2_pos)
        UMI2.collide <- .UMIcollide(UMIseq = UMI_seq2)
        
        names(UMI2.collide) <- names(chop_data$chopread)
        
        UMI2 <- DNAStringSet(as.character(UMI2.collide))
    }
    
    if(is.null(UMI_adaptor1_pos)){
        return(UMI2 = UMI)
    }else if(is.null(UMI_adaptor2_pos)){
        return(UMI1 = UMI)
    }else{
        return(list(UMI1 = UMI1, UMI2 = UMI2))
    }
        
}
    
    .UMIindelextract <- function(alignment, UMI_pos)
    {
        UMI_seq <- vector("list", length(alignment))
        for (i in 1:length(alignment)){
            shift <- sum(width(indel(subject(alignment[[i]])))[[1]][start(indel(subject(alignment[[i]])))[[1]] <= UMI_pos[2]])
            UMI_start <- UMI_pos[1] + shift 
            UMI_end <- UMI_pos[2] + shift
            UMI_seq[[i]] <- DNAStringSet(str_sub(pattern(alignment[[i]]),UMI_start, UMI_end))
        }
        
        return(UMI_seq)
    }
    
    .UMIcollide <- function(UMIseq)
    {
        UMI.collide <- vector("list", length(UMIseq))
        for (i in 1:length(UMI.collide)){
            cuts <- grep("[A-Z]", strsplit(as.character(UMIseq[[i]]), NULL)[[1]])
            x <- vector("list", length(cuts))
            for (a in 1:length(cuts)){
                x[a] <- strsplit(as.character(UMIseq[[i]]), NULL)[[1]][cuts[a]]
            }
            UMI.collide[i] <- paste(x, collapse = "")
        } 
        
        return(UMI.collide)
    }


