#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam bpstart bpstop bpisup bpmapply
qualityAlign <- function(sequences, reference, gapOpening=5, gapExtension=1, BPPARAM=SerialParam())
{
    by.core <- .parallelize(sequences, BPPARAM)
    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

	ref <- toupper(as.character(reference))
    out <- bpmapply(FUN=.align_general_internal, sequences=by.core, 
        MoreArgs=list(reference=ref, gap.opening=gapOpening, gap.extension=gapExtension),
        BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)

    output <- mapply(out, FUN=c)
    names(output) <- c("score", "reference", "query")
    output <- do.call(DataFrame, output)
    metadata(output) <- list(gapOpening=gapOpening, gapExtension=gapExtension, reference=reference)
    output
}

#' @importFrom Biostrings quality
.align_general_internal <- function(sequences, reference, gap.opening, gap.extension) {
    quals <- quality(sequences)
    .Call(cxx_general_align, sequences, quals, .create_encoding_vector(quals), 
            gap.opening, gap.extension, reference)
}
