#' @export
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom methods is
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
debarcodeReads <- function(align.stats, barcodes, position)
# Pulls out the barcodes, aligns them against all possible options,
# and reports the results.
#
# written by Aaron Lun
# created 19 September 2018
{
    go <- metadata(align.stats)$gapOpening 
    ge <- metadata(align.stats)$gapExtension
    ma <- metadata(align.stats)$match
    mm <- metadata(align.stats)$mismatch

    extracted <- umiExtract(align.stats, position=position)
    has.quality <- is(extracted, "QualityScaledDNAStringSet")
    all.args <- .setup_alignment_args(has.quality, go, ge, ma, mm, type="global")

    current.score <- rep(-Inf, length(barcodes))
    current.id <- rep(NA_integer_, length(barcodes))
    for (b in seq_along(barcodes)) {
        current <- .assign_qualities(barcodes[b])
        scores <- do.call(pairwiseAlignment, c(list(pattern=extracted, subject=current, scoreOnly=TRUE), all.args))
        keep <- scores > current.score
        current.id[keep] <- b
        current.score[keep] <- scores[keep]
    }

    DataFrame(barcode=current.id, score=current.score)
}
