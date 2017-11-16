
chopAndFilterReads <- function(align_data, names = NULL, essential1 = TRUE, essential2 = TRUE, align_score1 = NULL, align_score2 = NULL, align_pos1 = NULL, align_pos2 = NULL)
{
    if(essential1){
        id1 <- vector("logical", length(align_data$adaptor1))
        for (i in 1:length(align_data$adaptor1)){
            id1[i] <- start(pattern(align_data$adaptor1[[i]])) < align_pos1 & score(align_data$adaptor1[[i]]) >= align_score1
        }
    }else{
        id1 <- rep(TRUE, length(align_data$adaptor1))
    }
        
    if(essential2){
        id2 <- vector("logical", length(align_data$adaptor2))
        for (i in 1:length(align_data$adaptor2)){
            id2[i] <- nchar(align_data$reads[[i]])-end(pattern(align_data$adaptor2[[i]])) < align_pos2 & score(align_data$adaptor2[[i]]) >= align_score2
        }
    }else{
        id2 <- rep(TRUE, length(align_data$adaptor1))
    }

    id <- id1&id2
    reads <- align_data$reads[id]
    quality <- align_data$quality[id]

    if(essential1){
        adaptor1_end <- vector("integer", length = sum(id))
        chop_read1 <- vector("list", length = sum(id))
        align_adaptor1_filt <- vector("list", length = sum(id))
        for (i in which(id)){
            adaptor1_end[i] <- end(pattern(align_data$adaptor1[[i]]))
            align_adaptor1_filt[[i]] <- align_data$adaptor1[[i]]
        }
        chop_read1 <- subseq(reads, start = adaptor1_end+1)
        chop_qual1 <- subseq(quality, start = adaptor1_end+1)
        
    }else{
        adaptor1_end <- vector("integer", length = sum(id))
        chop_read1 <- reads
        chop_qual1 <- quality
        align_adaptor1_filt <- NULL
        
    }
    
    if(essential2){
        adaptor2_start <- vector("integer", length = sum(id))
        chop_read2 <- vector("list", length = sum(id))
        align_adaptor2_filt <- vector("list", length = sum(id))
        for (i in which(id)){
            adaptor2_start[i] <- start(pattern(align_data$adaptor2[[i]]))
            align_adaptor2_filt[[i]] <- align_data$adaptor2[[i]]
        }
        chop_read2 <- subseq(chop_read1, end = ((nchar(reads)-(102-adaptor2_start))-adaptor1_end))
        chop_qual2 <- subseq(chop_qual1, end = ((nchar(reads)-(102-adaptor2_start))-adaptor1_end))
        
        aligned_adaptors <- c(align_adaptor1_filt, align_adaptor2_filt)
        
    }else{
        chop_read2 <- chop_read1
        chop_qual2 <- chop_qual1
        aligned_adaptors <- align_adaptor1_filt
    }
    
    names(chop_read2) <- names[id]
    names(chop_qual2) <- names[id]
    #Read quality - if provided:
    
    return(list(chop_read2 = chopread, aligned_adaptors = adaptors, chop_qual2 = chopquality))
}
