\name{qualityAlign}
\alias{qualityAlign}

\title{Perform quality-based alignment}
\description{Align a quality-based sequence to a known reference.}

\usage{
qualityAlign(sequences, reference, gapOpening=5, gapExtension=1, 
    BPPARAM=SerialParam())
}

\arguments{
\item{sequences}{A \linkS4class{QualityScaledDNAStringSet} of observed barcode sequences, usually one of the subsequences from \code{\link{adaptorAlign}}.}
\item{reference}{A character vector or \linkS4class{DNAStringSet} containing the reference sequence.}
\item{gapOpening}{A numeric scalar specifying the gap opening penalty.}
\item{gapExtension}{A numeric scalar specifying the gap extension penalty.}
\item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying whether alignment should be parallelized.}
}

\details{
This provides an alternative to quality-aware alignment with \code{\link{pairwiseAlignment}} when one of the sequences is known. 
}

\value{
A \linkS4class{DataFrame} containing:
\describe{
\item{\code{score}:}{Numeric, the alignment score to the assigned barcode.}
\item{\code{reference}:}{Character, the alignment string for the reference sequence.}
\item{\code{sequence}:}{Character, the alignment string for each query sequence.}
}
Each row corresponds to an entry in \code{sequences}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{barcodeAlign}} and \code{\link{adaptorAlign}}, which operate on similar principles.
}

\examples{
ref <- "CAGCTGGACGATCAGCTCGATGCGACTGTGTGCATCG"
queries <- c("CAGGACGATCAGCTCGACGACTAGTGTGCATCGA",
    "CCTAAAGGACGATCAGCTCGATCTGTGTGCATCAAGA")

# Adding some qualities.
qquals <- lapply(nchar(queries), runif)
qquals <- do.call(c, lapply(qquals, as, "PhredQuality"))
queries <- QualityScaledDNAStringSet(queries, qquals)

# Assigning reads to individual barcodes.
align <- qualityAlign(queries, ref)
}
