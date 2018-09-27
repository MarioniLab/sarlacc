#' @export
#' @importFrom S4Vectors metadata
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom BiocParallel SerialParam
#' @importFrom methods is
getBarcodeThresholds <- function(aligned, error=0.01, BPPARAM=SerialParam())
# Computes thresholds for the barcode alignments,
# by comparison to scrambled sequences from the barcode segment of the read.
# 
# written by Aaron Lun
# created 27 September 2018
{
    seq <- align$sequence
    scrambled <- .scramble_input(seq, is(seq, "QualityScaledDNAStringSet"))
    bscores <- .align_to_barcodes(extracted, barcodes, metadata(align.stats), BPPARAM)
    .compute_threshold(aligned$score, bscores$score, error)
}
