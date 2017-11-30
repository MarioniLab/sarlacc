consensusReadSeq <- function(alignments, pseudo.count=1, min.coverage=0.6)
#ConsensusSequence by majority decision
#Returns List of consensus sequences for each cluster id and a list with the corresponding phredscores (which is named by the coverage of reads per consensus)
{
    Nalign <- nrow(alignments)
    consensus <- character(Nalign)
    phred <- list(Nalign)
    nucleotides <- c("A", "T", "C", "G")
    
    for (i in seq_along(alignments)) {
        current <- alignments[[i]]
        if (length(current)==1L) { 
            consensus[[i]] <- as.character(current)[1]
            phred[[i]] <- PhredQuality(rep(1/(1+pseudo.count), nchar(aln_char)))
        }
  
        # Keeping only the base positions with sufficient coverage. 
        x <- consensusMatrix(current) 
        x <- x[nucleotides,,drop=FALSE] # Why not x <- x[rownames(x)%in% nucleotides,]
        base.exists <- colSums(x) >= min.coverage * length(current)
        x <- x[,base.exists,drop=FALSE]

        # Defining the consensus sequence based on relative majority coverage
        # (at least one base must be above 25%, as we only keep 'x' in "ATCG")
        x <- t(x)
        majority.base <- max.col(x)
        conseq <- nucleotides[majority.base]
        consensus[i] <- paste(conseq, collapse="")
    
        # Computing the Phred score.
        n.chosen <- numeric(length(conseq))
        for (j in nucleotides) { 
            chosen <- conseq==j
            n.chosen[chosen] <- x[chosen,j]
        }
        phred[[i]] <- PhredQuality(1 - n.chosen/(rowSums(x)+pseudo.count))
    }

    return(list(sequence=DNAStringSet(consensus),
                quality=do.call(c, phred)))
}


# Majority decision does not work:
# [355,] 0 1 1 1
# [356,] 0 0 0 3
# [357,] 0 0 3 0
# [358,] 0 1 0 2
# [359,] 0 1 0 2
# [360,] 1 1 1 0
# -> 2 4 3 4 4 2