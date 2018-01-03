umiExtract <- function(align.stats, position=NULL) 
# Pull out UMIs from the adaptor alignment data.
#
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 27 November 2017    
{
    if (is.null(position)) { 
        curseq <- gsub("-", "", align.stats$adaptor[1])
        all.Ns <- gregexpr("N+", curseq)[[1]]
        position <- c(all.Ns, all.Ns+attr(all.Ns, "match.length")-1L)
    } else {
        if (length(position)!=2L || position[1] > position[2] || position[1] <= 0L) {
            stop("invalid 'position' vector")
        }
    }
    
    # Identifying the number of deletions in the adaptor before the start of the UMI.
    # This is the number of bases we have to shift to the 3' in the read alignment string. 
    read.bump <- .compute_position_bump(align.stats$adaptor, position)
    bumped.start <- read.bump$start + position[1]
    bumped.end <- read.bump$end + position[2]

    # Pulling out the UMI and cleaning it up.
    umi <- substr(align.stats$read, start=bumped.start, stop=bumped.end)
    umi <- gsub("-", "", umi)
    names(umi) <- rownames(align.stats)
    out <- DNAStringSet(umi)

    if (!is.null(align.stats$quality)) { 
        # Pulling out the Phred scores (this time, we need to know the deletions in the read,
        # as the Phred scores do not contain the deletions in the alignment string).
        read.unbump <- .compute_position_bump(align.stats$read, c(bumped.start, bumped.end))
        unbumped.start <- bumped.start - read.unbump$start[1]
        unbumped.end <- bumped.end - read.unbump$end[1]

        umi.qual <- subseq(align.stats$quality, start=unbumped.start, end=unbumped.end)
        out <- QualityScaledDNAStringSet(out, umi.qual)
    }
    return(out)
}

.compute_position_bump <- function(align.str, position) {
    all.dels <- gregexpr("-", align.str)
    bump.start <- bump.end <- integer(length(align.str))
    for (i in seq_along(bump.start)) { 
        dels <- as.integer(all.dels[[i]])
        if (length(dels)==1L && dels==-1L) { 
            next
        }
        true.pos <- dels - seq_along(dels)
        bump.start[i] <- sum(true.pos < position[1])
        bump.end[i] <- sum(true.pos < position[2])
    }
    return(list(start=bump.start, end=bump.end))
}

