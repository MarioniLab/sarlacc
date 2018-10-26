\name{expectedDist}
\alias{expectedDist}

\title{Calculate expected distances}
\description{Calculate expected distances between subsequences of the adaptor that should be identical across reads.}

\usage{
expectedDist(sequences, max.err=NA)
}

\arguments{
\item{sequences}{A \linkS4class{QualityScaledDNAStringSet} of read subsequences corresponding to constant regions of the adaptor.}
\item{max.err}{A numeric scalar specifying the maximum error probability above which bases will be masked.}
}

\details{
The aim is to provide an expectation for the distance for identical subsequences, given that all reads should originate from molecules with the same adaptor.
In this manner, we can obtain an appropriate threshold for \code{\link{umiGroup}} that accounts for sequencing and amplification errors.
We suggest extracting a subsequence from the interval next to the UMI region.
This ensures that the error rate in the extracted subsequence is as similar as possible to the UMI at that position on the read.

Pairwise Levenshtein distances are computed between all extracted sequences.
This is quite computationally expensive, so we only process a random subset of these sequences by setting \code{number}.
If \code{align.stats} contains quality scores, bases with error probabilities above \code{max.qual} are replaced with \code{N}s.
Any \code{N}s are treated as missing and will contribute a mismatch score of 0.5, even for matches to other \code{N}s.
}

\value{
A numeric vector of pairwise distances between sequences that \emph{should} be identical. 
}

\author{
Florian Bieberich,
with modifications by Aaron Lun
}

\seealso{
\code{\link{extractSubseq}} to extract a subsequence.
}

\examples{
constants <- c("ACTAGGAGA",
               "ACTACGACCA",
               "ACTACGATA",
               "ACACGACA")
expectedDist(constants)
}
