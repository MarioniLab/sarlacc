#' @export
#' @importFrom BiocParallel bplapply SerialParam
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom Biostrings quality encoding
#' @importFrom methods is new as
#' @importFrom S4Vectors metadata
#' @importClassesFrom S4Vectors DataFrame List
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
        by.groups <- list(seq_len(Nreads))
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
    do.mask <- !is.na(max.error)

    # Converting everything to a string, which is faster to work with.
    all.masked <- NULL
    all.reads <- as.character(reads)
    if (has.quals && !is.na(max.error)) {
        all.masked <- qualityMask(reads, max.error)
    } 

    # Pruning out singleton reads to avoid looping over them. 
    all.results <- vector("list", length(by.group))
    solo <- lengths(by.group)==1L
    if (any(solo)) {
        if (do.mask && keep.masked) {
            all.results[solo] <- all.masked[unlist(by.group[solo])]
        } else {
            all.results[solo] <- all.reads[unlist(by.group[solo])]
        }
    }

    # Running the multiple sequence alignment across multiple cores.
    multi <- !solo
    if (any(multi)) { 
        all.results[multi] <- bplapply(by.group[multi], FUN=.internal_msa, reads=all.reads, masked=all.masked, 
            ..., keep.masked=keep.masked, BPPARAM=BPPARAM)
    }

    # Collating the results into a single DF.
    out <- new("DataFrame", nrows=length(by.group))
    out$alignments <- as(all.results, "List")
    if (has.quals) {
        all.qual <- quality(reads)
        out$qualities <- as(lapply(by.group, function(idx) all.qual[idx]), "List")
    }
    rownames(out) <- names(by.group)
    return(out)
}

#' @importFrom Biostrings unmasked DNAStringSet
#' @importFrom muscle muscle
.internal_msa <- function(indices, reads, masked, keep.masked, ...) 
# Optimized MUSCLE call for speed.
{
    do.mask <- !is.null(masked)
    to.use <- if (do.mask) masked[indices] else reads[indices]
    cur.align <- muscle(DNAStringSet(to.use), ..., quiet=TRUE, diags=TRUE, maxiters=2)
    cur.align <- unmasked(cur.align)
    
    if (do.mask && !keep.masked) {
        cur.align <- .Call(cxx_unmask_alignment, cur.align, reads[indices])
    } else {
        cur.align <- as.character(cur.align)
    }
    cur.align
}
