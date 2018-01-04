multiReadAlign <- function(reads, groups, flip=NULL, min.qual=10, ...)
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
    if (!is.null(flip) && length(flip)!=Nreads) {
        stop("length of 'flip' and 'reads' should be the same")
    }

    by.group <- split(seq_len(Nreads), groups)
    msalign <- vector("list", length(by.group))
    names(msalign) <- names(by.group)
    reads <- .mask_bad_bases(reads, threshold=min.qual)

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

        cur.align <- muscle(cur.reads, ..., quiet=TRUE)
        msalign[[g]] <- DNAStringSet(cur.align)
    }
    return(msalign)
}

.mask_bad_bases <- function(incoming, threshold) {
    if (!is(incoming, "QualityScaledDNAStringSet")) {
        return(incoming)
    } else if (is.na(threshold)) { 
        return(as(incoming, "DNAStringSet"))        
    } else {
        enc <- encoding(quality(incoming))
        lowerbound <- names(enc)[min(which(enc>=threshold))]
        out <- .Call(cxx_mask_bad_bases, as.character(incoming),
                     as.character(quality(incoming)), 
                     max(width(incoming)), lowerbound)
        return(DNAStringSet(out))
    }
}
