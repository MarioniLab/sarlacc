multiReadAlign <- function(reads, groups, flip=NULL, min.qual=10, keep.masked=FALSE, ...)
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

    # Setting up some other constants.
    if (!is.null(flip) && length(flip)!=Nreads) {
        stop("length of 'flip' and 'reads' should be the same")
    }
    has.quals <- is(reads, "QualityScaledDNAStringSet") 
    do.mask <- !is.na(min.qual)

    by.group <- split(seq_len(Nreads), groups)
    msalign <- vector("list", length(by.group))
    names(msalign) <- names(by.group)

    for (g in names(by.group)) { 
        cur.reads <- reads[by.group[[g]]]
        if (length(cur.reads)==1L) {
            msalign[[g]] <- cur.reads
            next
        }

        if (!is.null(flip)) {
            curflip <- flip[by.group[[g]]]
            cur.reads[curflip] <- reverseComplement(cur.reads[curflip])
        }

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
            cur.align <- .unmask_bases(cur.align, to.use) 
        }

        # Storing the quality strings in the elementwise-metadata.
        if (has.quals) {
            elementMetadata(cur.align)$quality <- quality(cur.reads)
        }
        msalign[[g]] <- cur.align
    }
    return(msalign)
}

.mask_bad_bases <- function(incoming, threshold) {
    enc <- encoding(quality(incoming))
    lowerbound <- names(enc)[min(which(enc>=threshold))]
    out <- .Call(cxx_mask_bad_bases, incoming, quality(incoming), lowerbound)
    return(DNAStringSet(out))
}

.unmask_bases <- function(alignments, originals) {
    out <- .Call(cxx_unmask_bases, alignments, originals)
    return(DNAStringSet(out))
}
