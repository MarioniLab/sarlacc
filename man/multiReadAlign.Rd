\name{multiReadAlign}
\alias{multiReadAlign}

\title{Multiple read alignment}
\description{Perform a multiple sequence alignment of reads in the same UMI group, i.e., originating from the same DNA molecule.}

\usage{
multiReadAlign(reads, groups, max.error=0.1, keep.masked=FALSE, ..., 
    BPPARAM=SerialParam())
}

\arguments{
    \item{reads}{A DNAStringSet object containing read sequences.}
    \item{groups}{A factor specifying which reads belong to the same UMI group, e.g., from \code{\link{umiGroup}}.
    Alternatively, a list of indices specifying reads in the same group.}
    \item{max.error}{A numeric scalar specifying the maximum base-calling error,above which bases will be masked during MUSCLE.}
    \item{keep.masked}{A logical scalar indicating whether masked bases should remain masked in the output.}
    \item{...}{Additional arguments to pass to \code{\link{muscle}}.}
    \item{BPPARAM}{A BiocParallelParam object controlling how  paralellization should be performed.}
}

\details{
This function performs a multiple sequence alignment of reads from the same DNA molecule.
Reads should be grouped beforehand based on similarity in the UMI sequences, as described in \code{\link{umiGroup}}.
The aim is to use the alignments to form a consensus sequence using \code{\link{consensusReadSeq}}.
This allows us to obtain accurate sequences of the original molecule, overcoming the high error rate of ONT.

The MUSCLE aligner is used to perform efficient multiple sequence alignment of many reads.
Note that, for DNA sequences, MUSCLE assumes that all reads are on the same strand.
This should automatically be the case if \code{\link{adaptorAlign}} was used with two different adaptors on either end of the read.
However, if the same adaptors were used on both ends, the strand of the sequence is unknown and may require other information to discern.

MUSCLE is not aware of quality strings, so we mask low-quality bases prior to the multiple sequence alignment.
All bases with error probabilities above \code{max.error} are set to \code{"N"}.
This ensures that poor base calls do not compromise the alignment quality.
If \code{min.qual=NA}, no masking is performed.
All masks are removed in the output strings unless \code{keep.mask=TRUE}, in which case the strings may contain \code{"N"}'s.
}

\value{
A DataFrame containing \code{alignments}, a list of character vectors.
Each vector corresponds to an element of \code{groups} and contains the alignment strings for that group.

If \code{reads} is a QualityScaledDNAStringSet, the output DataFrame will also contain \code{qualities}, a list of lists of numeric vectors of base-calling error probabilities.
Each internal list corresponds to a MSA in \code{alignments}, and each vector corresponds to a sequence in the MSA.
}

\author{
Florian Bieberich,
with modifications from Aaron Lun
}

\seealso{
    \code{\link{umiGroup}},
    \code{\link{consensusReadSeq}},
    \code{\link{muscle}}
}

\references{
Edgar RC (2004). 
MUSCLE: multiple sequence alignment with high accuracy and high throughput. 
\emph{Nucleic Acids Res.}, 32, 5:1792-7.
}

\examples{
reads <- DNAStringSet(c("ACACTGGTTCAGGT", 
    "ACACGGTTCAGGT",
    "CGGACTGACACGGT",
    "CGGGCTGACACGGT"))

multiReadAlign(reads, c(1,1,2,2))

qreads <- QualityScaledDNAStringSet(reads,
    PhredQuality(c("23849723948733", 
    "*(&^&23498342",
    ")(*!@!343928&^",
    "2($*&$*A323128")))

multiReadAlign(qreads, c(1,1,2,2))
}


