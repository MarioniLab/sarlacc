chopAndFilterReads <- function(aligned, essential1 = TRUE, essential2 = TRUE, score1 = NULL, score2 = NULL) 
# This filters out reads that don't have essential adaptors aligning on either or both ends.
# We also chop out the adaptor sequences for furture use.    
{
    # Dropping reads if adaptor is essential but doesn't get detected.
    if(essential1){
        id1 <- aligned$adaptor1$score >= score1
    } else{
        id1 <- rep(TRUE, nrow(aligned$adaptor1))
    }

    if (essential2){
        id2 <- aligned$adaptor2$score >= score2
    }else{
        id2 <- rep(TRUE, nrow(align_data$adaptor1))
    }
    
    keep <- id1 & id2
    aligned$reads <- aligned$reads[keep]
    aligned$quality <- aligned$quality[keep]
    aligned$reversed <- aligned$reversed[keep]
    aligned$adaptor1 <- aligned$adaptor1[keep,]   
    aligned$adaptor2 <- aligned$adaptor2[keep,]   

    # Finding the cut points of each adaptor.
    # Score is filtered again to also mark adaptors that are not essential for chopping.
    start_point <- rep(1L, nrow(aligned$adaptor1))
    has1 <- aligned$adaptor1$score >= score1
    start_point[has1] <- aligned$adaptor1$end.pattern[has1] + 1L

    end_point <- width(aligned$reads)
    has2 <- aligned$adaptor2$score >= score2
    end_point[has2] <- aligned$adaptor2$start.pattern[has2] - 1L

    aligned$reads <- subseq(aligned$reads, start=start_point, end=end_point)
    aligned$quality <- subseq(aligned$quality, start=start_point, end=end_point)

    # Destroying pattern positional information, as this is no longer valid after chopping.
    if (essential1) { 
        aligned$adaptor1$start.pattern <- NULL    
        aligned$adaptor1$end.pattern <- NULL    
    } else {
        aligned$adaptor1 <- NULL
    }
    if (essential2) {
        aligned$adaptor2$start.pattern <- NULL    
        aligned$adaptor2$end.pattern <- NULL 
    } else {
        aligned$adaptor2 <- NULL
    }

    return(aligned)
}
