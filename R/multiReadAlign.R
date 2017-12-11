multiReadAlign <- function(reads, groups, flip=NULL, ...)
# Returns a DNAStringSet of multiple sequence alignments for the sequences from each cluster id
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
    by.group$"0" <- NULL
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

        cur.align <- muscle(cur.reads, ..., quiet=TRUE)
        msalign[[g]] <- DNAStringSet(cur.align)
    }
    return(msalign)
}
