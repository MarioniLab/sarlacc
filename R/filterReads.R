#' @export
filterReads <- function(aligned, score1, score2, essential1 = TRUE, essential2 = TRUE)
# This filters out reads that don't have essential adaptors aligning on either or both ends.
# We also chop out the adaptor sequences for future use.
#
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 24 November 2017
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
        id2 <- rep(TRUE, nrow(aligned$adaptor1))
    }

    keep <- id1 & id2
    aligned <- aligned[keep,]

    # Finding the cut points of each adaptor on the read sequence.
    # Score is filtered again to also mark adaptors that are not essential for chopping.
    # We also discard reads where adaptor 1 ends before adaptor 2 begins.
    start_point <- rep(1L, nrow(aligned$adaptor1))
    has1 <- aligned$adaptor1$score >= score1
    start_point[has1] <- aligned$adaptor1$end[has1] + 1L

    end_point <- aligned$read.width
    has2 <- aligned$adaptor2$score >= score2
    end_point[has2] <- aligned$adaptor2$end[has2] - 1L

    keep <- start_point < end_point
    aligned <- aligned[keep,]
    aligned$trimmed.start <- start_point[keep]
    aligned$trimmed.end <- end_point[keep]

    return(aligned)
}
