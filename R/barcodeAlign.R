#' @export
#' @importFrom S4Vectors metadata 
#' @importFrom BiocParallel SerialParam
barcodeAlign <- function(align.stats, barcodes, position, BPPARAM=SerialParam())
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
    all.args <- .setup_alignment_args(TRUE, param.list$gapOpening, param.list$gapExtension, type="global")
    all.args$scoreOnly <- TRUE
    by.core <- .parallelize(sequences, BPPARAM)

    current.score <- next.best <- rep(-Inf, length(sequences))
    current.id <- rep(NA_integer_, length(sequences))
    for (b in seq_along(barcodes)) {
        current <- .assign_qualities(barcodes[b])
        scores <- bpmapply(FUN=pairwiseAlignment, pattern=by.core, 
            MoreArgs=c(list(subject=current), all.args),
            BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
        scores <- unlist(scores)

        # Need to update both the best and the next best on record.
        keep <- scores > current.score
        second.keep <- !keep & scores > next.best

        current.id[keep] <- b
        next.best[keep] <- current.score[keep]
        current.score[keep] <- scores[keep]
        next.best[second.keep] <- scores[second.keep]
    }

    list(barcode=current.id, score=current.score, gap=current.score - next.best)
}
