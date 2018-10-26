#' @export
expectedDist <- function(sequences, max.err=NA)
# This computes the expected distance between a sequence of 
# the same length as the UMI.
#
# written by Florian Bieberich
# with modifications by Aaron Lun
{
    combined <- qualityMask(sequences, max.err)   
    .Call(cxx_compute_lev_masked, combined)
}
