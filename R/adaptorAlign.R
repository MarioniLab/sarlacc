#' @export
#' @importFrom Biostrings reverseComplement 
#' @importFrom S4Vectors DataFrame metadata<-
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom methods is
#' @importFrom BiocParallel SerialParam
adaptorAlign <- function(adaptor1, adaptor2, reads, tolerance=250, gapOpening=5, gapExtension=1, match=1, mismatch=0, BPPARAM=SerialParam())
# This function aligns both adaptors to the read sequence with the specified parameters,
# and returns the alignments that best match the sequence (with reverse complementing if necessary).
#    
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 7 November 2017    
{
    has.quality <- is(reads, "QualityScaledDNAStringSet")
    adaptor1 <- .assign_qualities(adaptor1, has.quality)
    adaptor2 <- .assign_qualities(adaptor2, has.quality)
    reads <- .assign_qualities(reads, has.quality)

    # Getting the start and (rc'd) end of the read.
    reads.out <- .get_front_and_back(reads, tolerance)
    reads.start <- reads.out$front
    reads.end <- reads.out$back

    # Performing the alignment of each adaptor to the start/end of the read.
    all.args <- .setup_alignment_args(has.quality, gapOpening, gapExtension, match, mismatch)
    align_start <- do.call(.bplalign, c(list(adaptor=adaptor1, reads=reads.start, BPPARAM=BPPARAM), all.args))
    align_end <- do.call(.bplalign, c(list(adaptor=adaptor2, reads=reads.end, BPPARAM=BPPARAM), all.args))
    align_revcomp_start <- do.call(.bplalign, c(list(adaptor=adaptor1, reads=reads.end, BPPARAM=BPPARAM), all.args))
    align_revcomp_end <- do.call(.bplalign, c(list(adaptor=adaptor2, reads=reads.start, BPPARAM=BPPARAM), all.args))
    
    # Standardizing the strand.
    strand.out <- .resolve_strand(align_start$score, align_end$score, align_revcomp_start$score, align_revcomp_end$score)
    is_reverse <- strand.out$reversed

    align_start[is_reverse,] <- align_revcomp_start[is_reverse,]
    align_end[is_reverse,] <- align_revcomp_end[is_reverse,]
    reads[is_reverse] <- reverseComplement(reads[is_reverse])

    details <- list(tolerance=tolerance, gapOpening=gapOpening, gapExtension=gapExtension, match=match, mismatch=mismatch)
    metadata(align_start) <- c(list(sequence=adaptor1), details)
    metadata(align_end) <- c(list(sequence=adaptor2), details)

    # Adjusting the reverse coordinates for the read length.
    old.start <- align_end$start
    old.end <- align_end$end
    align_end$start <- width(reads) - old.start + 1L
    align_end$end <- width(reads) - old.end + 1L
       
    rownames(align_start) <- rownames(align_end) <- names(reads) 
    output <- DataFrame(reads=reads, adaptor1=I(align_start), adaptor2=I(align_end), reversed=is_reverse, row.names=names(reads))
    return(output)
}

############################################
########### Internal functions #############
############################################

#' @importFrom Biostrings DNAStringSet PhredQuality
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom methods is
.assign_qualities <- function(truth, add.quality=TRUE)
# pairwiseAlignment() needs both strings to have qualities for quality-based pairwise alignment. 
# Obviously, the true reference has no qualities so we assign it the highest score possible.
{
    has.qual <- is(truth, "QualityScaledDNAStringSet")
    if (is.character(truth)) {
        truth <- DNAStringSet(truth)
    } 
    if (add.quality && !has.qual){
        qual <- PhredQuality(rep(100L, length(truth)))
        truth <- QualityScaledDNAStringSet(truth, qual)
    }
	truth  
}

#' @importFrom Biostrings reverseComplement
#' @importFrom XVector subseq
.get_front_and_back <- function(reads, tolerance) 
# Extracting the front part of the read (on the forward strand) 
# and the end of the read (on the reverse strand).
{
    tolerance <- pmin(tolerance, width(reads))
    reads.start <- subseq(reads, start = 1, width = tolerance)
    reads.end <- subseq(reads, end = width(reads), width = tolerance)
    reads.end <- reverseComplement(reads.end)
    return(list(front=reads.start, back=reads.end))
}

#' @importFrom Biostrings nucleotideSubstitutionMatrix
.setup_alignment_args <- function(has.quality, gapOpening, gapExtension, match, mismatch, type="local-global") {
    all.args <- list(type=type, gapOpening=gapOpening, gapExtension=gapExtension)
    if (!has.quality) { 
        all.args$substitutionMatrix <- nucleotideSubstitutionMatrix(match=match, mismatch=mismatch)
    } else {
        all.args$fuzzyMatrix <- nucleotideSubstitutionMatrix() # support IUPAC ambiguous codes in quality alignments. 
    }
    return(all.args)
}

.resolve_strand <- function(start.score, end.score, rc.start.score, rc.end.score) 
# This determines whether the read orientation needs to be flipped so that 
# adaptor 1 is at the start and adaptor 2 is at the end (see ?tuneAlignments 
# for an explanation of what the 'final.score' means).
{
    fscore <- pmax(start.score, 0) + pmax(end.score, 0)
    rscore <- pmax(rc.start.score, 0) + pmax(rc.end.score, 0)
    is.reverse <- fscore < rscore
    final.score <- ifelse(is.reverse, rscore, fscore)
    return(list(reversed=is.reverse, scores=final.score))
}

#' @importFrom Biostrings pairwiseAlignment quality pattern subject aligned unaligned
#' @importClassesFrom Biostrings QualityScaledDNAStringSet
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#' @importFrom utils head
#' @importFrom methods is
#' @importFrom BiocGenerics score start end
#' @importFrom S4Vectors DataFrame
#' @importFrom XVector subseq
.bplalign <- function(adaptor, reads, ..., scoreOnly=FALSE, BPPARAM=SerialParam())
# Align reads in parallel.
{
    n.cores <- bpnworkers(BPPARAM)
    bounds <- seq(from=1L, to=length(reads), length.out=n.cores+1L)
    ids <- findInterval(seq_along(reads), head(bounds, -1))
    by.core <- split(reads, ids)

    out <- bplapply(by.core, FUN=pairwiseAlignment, subject=adaptor, ..., scoreOnly=scoreOnly, BPPARAM=BPPARAM)
    if (scoreOnly) {
        return(unlist(out, use.names=FALSE))
    }

    has.quality <- is(reads, "QualityScaledDNAStringSet")
    collected <- vector("list", n.cores)
    for (i in seq_along(out)) {
        alignments <- out[[i]]
        read0 <- pattern(alignments)
        adaptor0 <- subject(alignments)
        extended <- .Call(cxx_get_aligned_sequence, aligned(adaptor0), as.character(unaligned(adaptor0)), start(adaptor0), end(adaptor0), aligned(read0))

        output <- DataFrame(score=score(alignments), adaptor=extended[[1]], read=extended[[2]], start=start(read0), end=end(read0))
        if (has.quality) {
            output$quality <- subseq(quality(by.core[[i]]), start=output$start, end=output$end)
        }
        collected[[i]] <- output
    }

    do.call(rbind, collected)
}
