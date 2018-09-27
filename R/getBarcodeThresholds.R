#' @export
#' @importFrom S4Vectors metadata
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom BiocParallel SerialParam
#' @importFrom methods is
getBarcodeThresholds <- function(baligned, error=0.01, BPPARAM=SerialParam())
# Computes thresholds for the barcode alignments,
# by comparison to scrambled sequences from the barcode segment of the read.
# 
# written by Aaron Lun
# created 27 September 2018
{
    seq <- baligned$sequence
    scrambled <- .scramble_input(seq, is(seq, "QualityScaledDNAStringSet"))
    bscores <- .align_to_barcodes(seq, metadata(baligned)$barcodes, metadata(baligned), BPPARAM)

    on.score <- .compute_threshold(baligned$score, bscores$score, error)
    on.diff <- .compute_threshold(baligned$gap, bscores$gap, error)
    c(score=on.score, gap=on.diff)
}
