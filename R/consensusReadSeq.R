#' @export
#' @importFrom Biostrings QualityScaledDNAStringSet DNAStringSet
#' @importFrom BiocParallel SerialParam bplapply
consensusReadSeq <- function(alignments, pseudo.count=1, min.coverage=0.6, BPPARAM=SerialParam())
# Create a consensus sequence for each MRA.
# 
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 27 November 2017    
{
    collected <- bplapply(alignments, .internal_consensus, pseudo.count=pseudo.count, min.coverage=min.coverage, BPPARAM=BPPARAM)        
    consensus <- unlist(lapply(collected, "[[", i=1))
    phred <- unname(lapply(collected, "[[", i=2))
    return(QualityScaledDNAStringSet(DNAStringSet(consensus), do.call(c, phred)))
}

#' @importFrom Biostrings PhredQuality
#' @importFrom S4Vectors elementMetadata
#' @importFrom methods as
.internal_consensus <- function(curalign, pseudo.count, min.coverage) { 
    quals <- elementMetadata(curalign)$quality
    has.quals <- !is.null(quals)

    # Skipping if we've only got one read in the alignment.
    if (length(curalign)==1L) { 
        consensus <- as.character(curalign)[1]
        if (has.quals) { 
            phred <- quals
        } else {
            phred <- PhredQuality(rep(1/(1+pseudo.count), nchar(consensus[i])))
        }
    } else {
        # Creating a consensus sequence that may or may not be Phred-aware.
        if (has.quals) {
            probs <- as.list(as(quals, "NumericList"))
            out <- .Call(cxx_create_consensus_quality, curalign, probs, min.coverage)
        } else {
            out <- .Call(cxx_create_consensus_basic, curalign, min.coverage, pseudo.count)
        }

        consensus <- out[[1]]
        phred <- PhredQuality(out[[2]])
    }

    list(consensus, phred)
}

