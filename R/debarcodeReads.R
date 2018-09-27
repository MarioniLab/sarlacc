#' @export
#' @importFrom S4Vectors metadata 
#' @importFrom BiocParallel SerialParam
debarcodeReads <- function(align.stats, barcodes, position, BPPARAM=SerialParam())
# Pulls out the barcodes, aligns them against all possible options,
# and reports the results.
#
# written by Aaron Lun
# created 19 September 2018
{
    extracted <- umiExtract(align.stats, position=position)
    bout <- .align_to_barcodes(extracted, barcodes, metadata(align.stats), BPPARAM)
    output <- DataFrame(sequence=extracted, bout, row.names=rownames(align.stats))
    metadata(output) <- c(metadata(align.stats)[c("gapOpening", "gapExtension", "match", "mismatch")], list(barcodes=barcodes))
    output
}

#' @importFrom methods is
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom S4Vectors DataFrame metadata<-
.align_to_barcodes <- function(sequences, barcodes, param.list, BPPARAM) {
    has.quality <- is(sequences, "QualityScaledDNAStringSet")
    all.args <- .setup_alignment_args(has.quality, param.list$gapOpening, param.list$gapExtension, 
            param.list$match, param.list$mismatch, type="global")
    all.args$reads <- sequences
    all.args$scoreOnly <- TRUE
    all.args$BPPARAM <- BPPARAM

    current.score <- next.best <- rep(-Inf, length(sequences))
    current.id <- rep(NA_integer_, length(sequences))
    for (b in seq_along(barcodes)) {
        current <- .assign_qualities(barcodes[b])
        scores <- do.call(.bplalign, c(list(adaptor=current), all.args))

        # Need to update both the best and the next best on record.
        keep <- scores > current.score
        second.keep <- !keep & scores > next.best

        current.id[keep] <- b
        next.best[keep] <- current.score[keep]
        current.score[keep] <- scores[keep]
        next.best[second.keep] <- scores[second.keep]
    }

    list(barcode=current.id, score=current.score, alternative=next.best)
}
