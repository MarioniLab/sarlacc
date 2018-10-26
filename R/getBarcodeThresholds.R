#' @export
#' @importFrom IRanges median mad
getBarcodeThresholds <- function(baligned, nmads=3)
# Computes thresholds for the barcode alignments,
# by comparison to scrambled sequences from the barcode segment of the read.
# 
# written by Aaron Lun
# created 27 September 2018
{
    med.s <- median(baligned$score)
    mad.s <- mad(baligned$score, center=med.s)
    med.g <- median(baligned$gap)
    mad.g <- mad(baligned$gap, center=med.g)
    c(score=med.s - mad.s*nmads, gap=med.g - mad.g*nmads)
}
