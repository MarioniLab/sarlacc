multiSeqAlign <- function(reads, groups, ...)
# Returns a DNAStringSet of multiple sequence alignments for the sequences from each cluster id
{
    Nreads <- length(reads)
    if (missing(groups)) {
        groups <- rep(1L, Nreads)
    }
    if (length(groups)!=Nreads) {
        stop("length of 'reads' and 'groups' should be the same")
    }   

    by.group <- split(seq_len(Nreads), groups)
    msalign <- vector("list", length(by.group))
    names(msalign) <- names(by.group)

    for (g in names(by.group)) { 
        cur.reads <- reads[by.group[[g]]]
        cur.align <- muscle(cur.reads, ..., quiet=TRUE)
        msalign[[g]] <- DNAStringSet(cur.align)
    }
    return(msalign)
}
