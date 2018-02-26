\name{multiReadAlign}
\alias{multiReadAlign}

\title{Multiple read alignment}
\description{Perform a multiple sequence alignment of reads in the same UMI group, i.e., originating from the same DNA molecule.}

\usage{
multiReadAlign(reads, groups, min.qual=10, keep.masked=FALSE, ...)
}

\arguments{
    \item{reads}{A DNAStringSet object containing read sequences.}
    \item{groups}{A factor specifying which reads belong to the same UMI group, e.g., from \code{\link{umiGroup}}.}
    \item{min.qual}{An integer scalar specifying the minimum quality score, below which bases will be masked during MUSCLE.}
    \item{keep.masked}{A logical scalar indicating whether masked bases should remain masked in the output.}
    \item{...}{Additional arguments to pass to \code{\link{muscle}}.}
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
All bases with quality scores (i.e., the \code{-10*log10(error)} probabilities) below \code{min.qual} are set to \code{"N"}.
This ensures that poor base calls do not compromise the alignment quality.
If \code{min.qual=NA}, no masking is performed.
All masks are removed in the output strings unless \code{keep.mask=TRUE}, in which case the strings may contain \code{"N"}'s.
}

\value{
A list of DNAStringSet objects.
Each DNAStringSet contains alignment strings for all reads in one UMI group.
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

multiReadAlign(reads, c(1,1,2,2), 
    flip=c(TRUE, TRUE, FALSE, FALSE))
}


