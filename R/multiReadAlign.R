#' @export
#' @importFrom BiocParallel bplapply SerialParam
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom Biostrings quality encoding
#' @importFrom methods is
multiReadAlign <- function(reads, groups, min.qual=10, keep.masked=FALSE, ..., BPPARAM=SerialParam())
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
    do.mask <- !is.na(min.qual)
    if (has.quals) {
        if (is.na(min.qual)) { 
            enc <- encoding(quality(reads))
            lowerbound <- names(enc)[min(which(enc>=min.qual))] 
        } else {
            lowerbound <- NA_character_
        }
    }

    # Converting everything to a string, which is faster to work with.
    all.masked <- all.qual <- NULL
    all.reads <- as.character(reads)
    if (has.quals) {
        all.qual <- as.character(quality(reads))
        if (do.mask) {
            all.masked <- .Call(cxx_mask_bad_bases, all.reads, all.qual, min.enc)
        }
    } 

    # Running the multiple sequence alignment across multiple cores.
    multi.res <- bplapply(by.group, FUN=.internal_msa, reads=all.reads, qualities=all.qual, masked=all.masked, 
        ..., keep.masked=keep.masked, BPPARAM=BPPARAM)

    # Collating the results into a single DF.
    out <- new("DataFrame", nrows=length(by.group))
    out$alignments <- multi.res 
    if (has.quals) {
        out$qualities <- lapply(by.group, function(x) all.qual[idx]) 
    }
    rownames(out) <- names(by.group)
    return(out)
}

#' @importFrom Biostrings unmasked
#' @importFrom muscle muscle
.internal_msa <- function(indices, reads, qualities, masked, keep.masked, ...) {
    do.mask <- !is.null(masked)
    if (length(indices)==1L) { 
        if (keep.masked) {
            return(masked[indices])
        } else {
            return(reads[indices])
        }
    }

    to.use <- if (do.mask) masked[indices] else reads[indices]
    cur.align <- muscle(DNAStringSet(to.use), ..., quiet=TRUE)
    cur.align <- unmasked(cur.align)
    
    if (do.mask && !keep.masked) {
        cur.align <- .Call(cxx_unmask_bases, cur.align, reads[indices])
    } else {
        cur.align <- as.character(cur.align)
    }
    cur.align
}
