#' @export
#' @importFrom BiocParallel bplapply SerialParam
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom Biostrings quality encoding
#' @importFrom methods is new
#' @importFrom S4Vectors metadata
#' @importClassesFrom S4Vectors DataFrame
multiReadAlign <- function(reads, groups, max.error=NA, keep.masked=FALSE, ..., BPPARAM=SerialParam())
# Returns a DNAStringSet of multiple sequence alignments for the sequences from each cluster id.
# MUSCLE itself is not quality-aware, so we help it out by masking low-quality bases beforehand.
#
# written by Florian Bieberich
# with modifications from Aaron Lun
# created 27 November 2017
{
    Nreads <- length(reads)
    if (missing(groups)) {
        by.group <- list(seq_len(Nreads))
    } else if (is.list(groups)) {
        by.group <- groups
    } else {
        if (length(groups)!=Nreads) {
            stop("length of 'reads' and 'groups' should be the same")
        }
        by.group <- split(seq_len(Nreads), groups)
    }

    # Setting quality-related parameters.
    has.quals <- is(reads, "QualityScaledDNAStringSet")
    if (has.quals) {
        all.qual <- as.list(as(quality(reads), "NumericList"))
    } else {
        all.qual <- NULL
    }

    # Setting up and running the MSA.
    dna5 <- c("A", "C", "G", "T", "N")
    submat <- nucleotideSubstitutionMatrix()[dna5, dna5, drop=FALSE]
    gapOpening <- 5
    gapExtension <- 1
    all.results <- .internal_msa(reads, by.group, max.error, all.qual, submat, -gapOpening, -gapExtension)

    # Collating the results into a single DF.
    out <- new("DataFrame", nrows=length(by.group))
    out$alignments <- all.results 
    if (has.quals) {
        out$qualities <- lapply(by.group, function(idx) all.qual[idx]) 
    }
    rownames(out) <- names(by.group)
    return(out)
}

.internal_msa <- function(sequences, read.groups, max.err, qualities, submat, opening, extension) {
    .Call(cxx_quick_msa, sequences, read.groups, max.err, qualities, submat, opening, extension) 
}
