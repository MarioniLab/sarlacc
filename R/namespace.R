#' @import Biostrings
#' @import XVector
#' @import IRanges
#' @import methods
#' @import S4Vectors
#' @import GenomicRanges
#' @importFrom Rcpp sourceCpp
#' @useDynLib sarlacc, .registration=TRUE, .fixes="cxx_"
NULL

# Note: need XVector, IRanges and methods "Imports" for the Biostrings C API.
# Note: need Biostrings, S4Vectors and GenomicRanges "Depends", as these are returned to the user.
