# Takes beginning of chopped read and computes the levenshtein distance to the earlier extracted UMI sequence. Acts as a negative control to levExpect.
# Gives information about similarity between random sequences.

# reads <- chop_data$reads
# umi <- UMI1

negExpect <- function(umi, reads, x = 100)
{
    reads <- reads[1:x,]
    umi <- umi[1:x,]
    
    # length <- round(mean(width(UMI1)), digits = 0)
    
    a <- DNAStringSet(substr(reads, start=1, stop=width(UMI1)))
    b <- umi
    
    lev.dist <- as.vector(stringDist(c(a,b), method = "levenshtein"))
    
    return(c(mean(lev.dist),sd(lev.dist)))
    
}
