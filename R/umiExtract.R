#' @export
#' @importFrom Biostrings DNAStringSet QualityScaledDNAStringSet
#' @importFrom XVector subseq
umiExtract <- function(align.stats, position=NULL) 
# Pull out UMIs from the adaptor alignment data.
#
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 27 November 2017    
{
    if (nrow(align.stats)==0L){ 
        stop("'align.stats' should have non-zero rows")
    }
    position <- .guess_umi_position(position, align.stats$adaptor[1])
    
    # Identifying the number of deletions in the adaptor before the start of the UMI.
    # This is the number of bases we have to shift to the 3' in the read alignment string. 
    read.bump <- .Call(cxx_count_gaps_by_base, align.stats$adaptor, position[1], position[2])
    bumped.start <- position[1] + read.bump[[1]]
    bumped.end <- position[2] + read.bump[[2]]

    # Pulling out the UMI and cleaning it up.
    umi <- substr(align.stats$read, start=bumped.start, stop=bumped.end)
    umi <- gsub("-", "", umi)
    names(umi) <- rownames(align.stats)
    out <- DNAStringSet(umi)

    if (!is.null(align.stats$quality)) { 
        # Pulling out the Phred scores (this time, we need to know the deletions in the read,
        # as the Phred scores do not contain the deletions in the alignment string).
        read.unbump <- .Call(cxx_count_gaps_by_align, align.stats$read, bumped.start, bumped.end)
        unbumped.start <- bumped.start - read.unbump[[1]]
        unbumped.end <- bumped.end - read.unbump[[2]]
        umi.qual <- subseq(align.stats$quality, start=unbumped.start, stop=unbumped.end)
        out <- QualityScaledDNAStringSet(out, umi.qual)
    }
    return(out)
}

.guess_umi_position <- function(position, adaptor) { 
    if (is.null(position)) { 
        curseq <- gsub("-", "", adaptor)
        all.Ns <- gregexpr("N+", curseq)[[1]]
        position <- c(all.Ns, all.Ns+attr(all.Ns, "match.length")-1L)
    } else {
        if (length(position)!=2L || position[1] > position[2] || position[1] <= 0L) {
            stop("invalid 'position' vector")
        }
    }
    return(position)
}
