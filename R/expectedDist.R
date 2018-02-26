#' @export
#' @importFrom Biostrings DNAStringSet QualityScaledDNAStringSet
#' @importFrom methods is
expectedDist <- function(align.stats, position=NULL, offset=NULL, number=100, get.seq=FALSE, min.qual=10)
# This computes the expected distance between a sequence of 
# the same length as the UMI.     
#
# written by Florian Bieberich
# with modifications by Aaron Lun
{
    if (nrow(align.stats)==0L) {
        stop("'align.stats' should have non-zero rows")
    }
    position <- .guess_umi_position(position, align.stats$adaptor[1])
    align.stats <- head(align.stats, number)

    # Figuring out where we're extracting from.
    umi.len <- position[2] - position[1] + 1L
    if (is.null(offset)) {
        offset <- floor(umi.len/2L)
    } else if (offset < 0 || offset > umi.len) {
        stop("offset must be an integer in [0, UMI length]")        
    }
    front.half <- position[1] - c(offset, 1L)
    back.half <- position[2] + c(1L, umi.len - offset)

    # Doing the extraction, unless the 'offset' is defined to avoid the back or front.
    if (offset!=0L) { 
        front.seq <- umiExtract(align.stats, position=front.half)
    } else {
        front.seq <- NULL 
    }
    if (offset!=umi.len) {
        back.seq <- umiExtract(align.stats, position=back.half)
    } else {
        back.seq <- subseq(front.seq, 1, 0)
    }
    if (is.null(front.seq)) { # needs to be done after back.seq is defined.
        front.seq <- subseq(back.seq, 1, 0)
    }

    # Combining the front and back sequences together.
    combined <- DNAStringSet(paste0(front.seq, back.seq))
    if (is(front.seq, "QualityScaledDNAStringSet")) {
        combined.qual <- paste0(quality(front.seq), quality(back.seq))
        combined <- QualityScaledDNAStringSet(combined, combined.qual)
    }

    if (get.seq) {
        return(combined)
    } 
    combined <- .safe_masker(combined, threshold=min.qual)    
    return(.Call(cxx_compute_lev_masked, combined))
}
