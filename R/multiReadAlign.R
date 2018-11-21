#' @export
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom methods is new as
#' @importClassesFrom S4Vectors DataFrame List
#' @importFrom Biostrings quality
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
multiReadAlign <- function(reads, groups, max.error=NA, match=0, mismatch=-1, gapOpening=5, gapExtension=1, bandwidth=100, BPPARAM=SerialParam())
# Returns multiple sequence alignments for the sequences from each cluster id.
# The aligner itself is not quality-aware, so we help it out by masking low-quality 
# bases beforehand if 'max.error' is set to some non-'NA' value.
#
# written by Florian Bieberich
# with modifications from Aaron Lun
# created 27 November 2017
{
    Nreads <- length(reads)
    if (missing(groups)) {
        by.groups <- list(seq_len(Nreads))
    } else if (is.list(groups)) {
        by.group <- groups
    } else {
        if (length(groups)!=Nreads) {
            stop("length of 'reads' and 'groups' should be the same")
        }
        by.group <- split(seq_len(Nreads), groups)
    }

    # Running the multiple sequence alignment across multiple cores.
    by.core <- .parallelize(by.group, BPPARAM)
    all.results <- bplapply(by.core, FUN=.internal_msa, reads=reads, match=match, mismatch=mismatch, 
        gapOpening=gapOpening, gapExtension=gapExtension, bandwidth=bandwidth, BPPARAM=BPPARAM)
    all.results <- unlist(all.results, recursive=FALSE)

    # Collating the results into a single DF.
    out <- new("DataFrame", nrows=length(by.group))
    out$alignments <- as(all.results, "List")
    rownames(out) <- names(by.group)

    if (is(reads, "QualityScaledDNAStringSet")) {
        all.qual <- quality(reads)
        out$qualities <- as(lapply(by.group, function(idx) all.qual[idx]), "List")
    }
    return(out)
}

.internal_msa <- function(groups, reads, match, mismatch, gapOpening, gapExtension, bandwidth) {
    .Call(cxx_quick_msa, groups, reads, match, mismatch, -gapOpening, -gapExtension, bandwidth)
}
