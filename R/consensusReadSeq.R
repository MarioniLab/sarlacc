#' @export
#' @importFrom Biostrings QualityScaledDNAStringSet DNAStringSet
#' @importFrom BiocParallel SerialParam bpmapply
consensusReadSeq <- function(alignments, pseudo.count=1, min.coverage=0.6, BPPARAM=SerialParam())
# Create a consensus sequence for each MRA.
#
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 27 November 2017
{
    aln <- alignments$alignments
    qual <- alignments$qualities
    if (is.null(qual)) {
        qual <- vector("list", nrow(alignments))
    }

    # Handling singletons separately.
    solo <- lengths(aln)==1L
    all.seq <- all.qual <- vector("list", length(aln))
    if (any(solo)) {
        solo.aln <- aln[solo]
        all.seq[solo] <- solo.aln
        if (has.quals) {
            all.qual[solo] <- qual[solo]
        } else {
            constant <- rep(1/(1+pseudo.count))
            all.seq[solo] <- lapply(solo.aln, function(seq) { rep(constant, nchar(seq)) })
        }
    }

    # Handling the multiples.
    multi <- !solo
    if (any(multi)) {
        collected <- bpmapply(.internal_consensus, aln[multi], qual[multi], 
            MoreArgs=list(pseudo.count=pseudo.count, min.coverage=min.coverage), 
            SIMPLIFY=FALSE, BPPARAM=BPPARAM)
        all.seq[multi] <- unlist(lapply(collected, "[[", i=1))
        all.qual[multi] <- unname(lapply(collected, "[[", i=2))
    }

    all.seq <- unlist(all.seq)
    all.qual <- do.call(c, lapply(all.qual, PhredQuality))
    return(QualityScaledDNAStringSet(DNAStringSet(all.seq), all.qual))
}

#' @importFrom Biostrings PhredQuality
.internal_consensus <- function(alignment, qualities, pseudo.count, min.coverage) {
    has.quals <- !is.null(qualities)
    # Creating a consensus sequence that may or may not be Phred-aware.
    if (has.quals) {
        out <- .Call(cxx_create_consensus_quality, alignment, min.coverage, qualities)
    } else {
        out <- .Call(cxx_create_consensus_basic, alignment, min.coverage, pseudo.count)
    }
    return(out)
}

