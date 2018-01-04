\name{multiReadAlign}
\alias{multiReadAlign}

\title{Multiple read alignment}
\description{Perform a multiple sequence alignment of reads in the same UMI group, i.e., originating from the same DNA molecule.}

\usage{
multiReadAlign(reads, groups, flip=FALSE, threshold=10, ...)
}

\arguments{
    \item{reads}{A DNAStringSet object containing read sequences.}
    \item{groups}{A factor specifying which reads belong to the same UMI group, e.g., from \code{\link{umiGroup}}.}
    \item{flip}{A logical vector specifying whether each read should be reverse-complemented.}
    \item{threshold}{An integer scalar specifying the minimum quality score, below which bases will be masked during MUSCLE.}
    \item{...}{Additional arguments to pass to \code{\link{muscle}}.}
}

\details{
This function performs a multiple sequence alignment of reads from the same DNA molecule.
Reads should be grouped beforehand based on similarity in the UMI sequences, as described in \code{\link{umiGroup}}.
The aim is to use the alignments to form a consensus sequence using \code{\link{consensusReadSeq}}.
This allows us to obtain accurate sequences of the original molecule, overcoming the high error rate of ONT.

The MUSCLE aligner is used to perform efficient multiple sequence alignment of many reads.
However, for DNA sequences, MUSCLE assumes that all reads are on the same strand.
Reads on the other strand should be indicated as such via \code{flip} so that their reverse complement sequences are used in the multiple sequence alignment.
In general, this is only necessary if the same adaptors were used on both ends such that the strand of the sequence is unknown.

MUSCLE is not aware of quality strings, so we mask low-quality bases prior to the multiple sequence alignment.
All bases with quality scores (i.e., the \code{-10*log10()} error probabilities) below \code{threshold} are set to \code{"N"}.
This ensures that poor base calls do not compromise the alignment quality.
If \code{threshold=NA}, no masking is performed.
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


