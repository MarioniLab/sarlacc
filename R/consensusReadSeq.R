#' @export
#' @importFrom Biostrings QualityScaledDNAStringSet DNAStringSet
#' @importFrom BiocParallel SerialParam bpmapply
#' @importFrom Biostrings PhredQuality
consensusReadSeq <- function(alignments, pseudo.count=1, min.coverage=0.6, BPPARAM=SerialParam())
# Create a consensus sequence for each MRA.
#
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 27 November 2017
{
    aln <- alignments$alignments
    qual <- alignments$qualities

    has.quals <- !is.null(qual)
    if (has.quals) {
        error.probs <- .create_encoding_vector(qual[[1]])
        out <- .Call(cxx_create_consensus_quality_loop, aln, min.coverage, qual, error.probs)
    } else {
        out <- .Call(cxx_create_consensus_basic_loop, aln, min.coverage, pseudo.count)
    }

    return(QualityScaledDNAStringSet(DNAStringSet(out[[1]]), PhredQuality(out[[2]])))
}
