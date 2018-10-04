#' @export
#' @importFrom Biostrings encoding
#' @importFrom methods is 
#' @importClassesFrom Biostrings DNAStringSet QualityScaledDNAStringSet
qualityMask <- function(seq, max.err) {
    has.quals <- is(seq, "QualityScaledDNAStringSet")
    if (!is.na(max.err) && has.quals) {
        quals <- quality(seq)
        error.probs <- .create_encoding_vector(quals)
        seq <- .Call(cxx_mask_bad_bases, seq, quals, error.probs, max.err)
    } else if (has.quals) {
        seq <- as(seq, "DNAStringSet")
    }
    return(seq)
}

#' @importFrom Biostrings encoding
#' @importFrom methods as
.create_encoding_vector <- function(quals) {
    enc <- encoding(quals)

    all.scores <- paste(names(enc), collapse="")
    all.scores <- as(all.scores, class(quals))
    error.probs <- as.numeric(all.scores)
    names(error.probs) <- names(enc)

    error.probs
}
