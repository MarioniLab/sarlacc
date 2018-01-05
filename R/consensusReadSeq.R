consensusReadSeq <- function(alignments, pseudo.count=1, min.coverage=0.6)
# Create a consensus sequence for each MRA.
# 
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 27 November 2017    
{
    Nalign <- length(alignments)
    consensus <- character(Nalign)
    phred <- vector("list", Nalign)
    
    for (i in seq_along(alignments)) {
        current <- alignments[[i]]
        quals <- elementMetadata(current)$quality
        has.quals <- !is.null(quals)

        # Skipping if we've only got one read in the alignment.
        if (length(current)==1L) { 
            consensus[[i]] <- as.character(current)[1]
            if (has.quals) { 
                phred[[i]] <- quals
            } else {
                phred[[i]] <- PhredQuality(rep(1/(1+pseudo.count), nchar(aln_char)))
            }
            next
        }

        # Creating a consensus sequence that may or may not be Phred-aware.
        if (has.quals) {
            probs <- as.list(as(quals, "NumericList"))
            out <- .Call(cxx_create_consensus_quality, current, probs, min.coverage)
        } else {
            out <- .Call(cxx_create_consensus_basic, current, min.coverage, pseudo.count)
        }

        consensus[i] <- out[[1]]
        phred[[i]] <- PhredQuality(out[[2]])
    }

    return(QualityScaledDNAStringSet(DNAStringSet(consensus), do.call(c, phred)))
}

