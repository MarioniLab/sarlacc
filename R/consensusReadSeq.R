#' @export
#' @importFrom Biostrings QualityScaledDNAStringSet DNAStringSet
#' @importFrom BiocParallel SerialParam bpmapply
consensusReadSeq <- function(alignments, pseudo.count=1, min.coverage=0.6, BPPARAM=SerialParam())
# Create a consensus sequence for each MRA.
#
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 27 November 2017
{
    aln <- alignments$alignments
    qual <- alignments$qualities
    if (is.null(qual)) {
        qual <- vector("list", nrow(alignments))
    }

    collected <- bpmapply(.internal_consensus, aln, qual, 
        MoreArgs=list(pseudo.count=pseudo.count, min.coverage=min.coverage), 
        SIMPLIFY=FALSE, BPPARAM=BPPARAM)

    consensus <- unlist(lapply(collected, "[[", i=1))
    phred <- unname(lapply(collected, "[[", i=2))
    return(QualityScaledDNAStringSet(DNAStringSet(consensus), do.call(c, phred)))
}

#' @importFrom Biostrings PhredQuality
.internal_consensus <- function(alignment, qualities, pseudo.count, min.coverage) {
    has.quals <- !is.null(qualities)

    # Skipping if we've only got one read in the alignment.
    if (length(alignment)==1L) {
        if (has.quals) {
            phred <- PhredQuality(qualities[[1]])
        } else {
            if (nchar(alignment)) {
                phred <- PhredQuality(rep(1/(1+pseudo.count), nchar(alignment)))
            } else {
                phred <- PhredQuality("")
            }
        }
        return(list(alignment, phred))
    }

    # Creating a consensus sequence that may or may not be Phred-aware.
    if (has.quals) {
        out <- .Call(cxx_create_consensus_quality, alignment, min.coverage, qualities)
    } else {
        out <- .Call(cxx_create_consensus_basic, alignment, min.coverage, pseudo.count)
    }

    out[[2]] <- PhredQuality(out[[2]])
    return(out)
}

