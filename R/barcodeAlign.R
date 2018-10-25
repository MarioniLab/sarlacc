#' @export
#' @importFrom S4Vectors metadata metadata<- DataFrame
#' @importFrom BiocParallel SerialParam bpstart bpstop bpisup bpmapply
barcodeAlign <- function(sequences, barcodes, gapOpening=5, gapExtension=1, BPPARAM=SerialParam())
# Pulls out the barcodes, aligns them against all possible options,
# and reports the results.
#
# written by Aaron Lun
# created 19 September 2018
{
    by.core <- .parallelize(sequences, BPPARAM)
    current.score <- next.best <- rep(-Inf, length(sequences))
    current.id <- rep(NA_integer_, length(sequences))

    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    for (b in seq_along(barcodes)) {
        current <- toupper(as.character(barcodes[b]))
        scores <- bpmapply(FUN=.align_BA_internal, barcodes=by.core, 
            MoreArgs=list(reference=barcodes[b], gap.opening=gapOpening, gap.extension=gapExtension),
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

    output <- DataFrame(barcode=current.id, score=current.score, gap=current.score - next.best)
    metadata(output) <- list(gapOpening=gapOpening, gapExtension=gapExtension, barcodes=barcodes)
    output
}

#' @importFrom Biostrings quality
.align_BA_internal <- function(barcodes, reference, gap.opening, gap.extension) {
    quals <- quality(barcodes)
    .Call(cxx_barcode_align, barcodes, quals, .create_encoding_vector(quals), 
            gap.opening, gap.extension, reference)
}
