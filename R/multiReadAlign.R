#' @export
#' @importFrom BiocParallel bplapply SerialParam
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
        groups <- rep(1L, Nreads)
    }
    if (length(groups)!=Nreads) {
        stop("length of 'reads' and 'groups' should be the same")
    }

    by.group <- split(seq_len(Nreads), groups)
    msalign <- bplapply(by.group, FUN=.internal_msa, reads=reads, min.qual=min.qual,
        keep.masked=keep.masked, ..., BPPARAM=BPPARAM)

    return(msalign)
}

#' @importClassesFrom Biostrings DNAStringSet QualityScaledDNAStringSet
#' @importFrom Biostrings unmasked quality
#' @importFrom S4Vectors elementMetadata<-
#' @importFrom methods is as
#' @importFrom muscle muscle
.internal_msa <- function(reads, group, min.qual, keep.masked, ...) {
    cur.reads <- reads[group]
    has.quals <- is(cur.reads, "QualityScaledDNAStringSet")
    do.mask <- !is.na(min.qual)

    if (length(cur.reads)==1L) {
        if (has.quals) {
            if (do.mask && keep.masked) {
                cur.align <- .mask_bad_bases(cur.reads, threshold=min.qual)
            } else {
                cur.align <- as(cur.reads, "DNAStringSet")
            }
        } else {
            cur.align <- cur.reads
        }
    } else {
        # Performing the MSA on potentially masked reads.
        if (has.quals) {
            if (do.mask) {
                to.use <- .mask_bad_bases(cur.reads, threshold=min.qual)
            } else {
                to.use <- as(cur.reads, "DNAStringSet")
            }
        } else {
            to.use <- cur.reads
        }
        cur.align <- muscle(to.use, ..., quiet=TRUE)
        cur.align <- unmasked(cur.align)

        # Unmasking the bases in the alignments.
        if (has.quals && do.mask && !keep.masked) {
            cur.align <- .unmask_bases(cur.align, cur.reads)
        }
    }

    # Storing the quality strings in the elementwise-metadata.
    if (has.quals) {
        elementMetadata(cur.align)$quality <- quality(cur.reads)
    }
    cur.align
}

#' @importFrom Biostrings encoding quality DNAStringSet
.mask_bad_bases <- function(incoming, threshold) {
    enc <- encoding(quality(incoming))
    lowerbound <- names(enc)[min(which(enc>=threshold))]
    out <- .Call(cxx_mask_bad_bases, incoming, quality(incoming), lowerbound)
    return(DNAStringSet(out))
}

#' @importFrom Biostrings DNAStringSet
.unmask_bases <- function(alignments, originals) {
    out <- .Call(cxx_unmask_bases, alignments, originals)
    return(DNAStringSet(out))
}

#' @importFrom Biostrings QualityScaledDNAStringSet DNAStringSet
#' @importFrom methods is
.safe_masker <- function(incoming, threshold) {
    if (is(incoming, "QualityScaledDNAStringSet") && !is.na(threshold)) {
        return(.mask_bad_bases(incoming, threshold))
    } else {
        return(as(incoming, "DNAStringSet"))
    }
}




