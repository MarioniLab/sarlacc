\name{multiReadAlign}
\alias{multiReadAlign}

\title{Multiple read alignment}
\description{Perform a multiple sequence alignment of reads in the same UMI group, i.e., originating from the same DNA molecule.}

\usage{
multiReadAlign(reads, groups, max.error=NA, match=0, mismatch=-1, 
    gapOpening=5, gapExtension=1, bandwidth=100, BPPARAM=SerialParam())
}

\arguments{
    \item{reads}{A \linkS4class{DNAStringSet} or \linkS4class{QualityScaledDNAStringSet} object containing read sequences.}
    \item{groups}{A factor specifying which reads belong to the same UMI group, e.g., from \code{\link{umiGroup}}.
    Alternatively, a list of indices specifying reads in the same group.}
    \item{max.error}{A numeric scalar specifying the maximum base-calling error, above which bases will be masked during alignment.
    Currently ignored.}
    \item{match, mismatch}{Numeric scalars specifying the match and mismatch scores.}
    \item{gapOpening, gapExtension}{Numeric scalars specifying the gap opening and extension penalties.}
    \item{bandwidth}{Integer scalar specifying the bandwidth for banded alignment, 
        corresponding to the maximum net difference in the numbers of insertions and deletions.}
    \item{BPPARAM}{A \linkS4class{BiocParallelParam} object controlling how  paralellization should be performed.}
}

\details{
This function performs a multiple sequence alignment of reads from the same DNA molecule.
Reads should be grouped beforehand based on similarity in the UMI sequences, as described in \code{\link{umiGroup}}.
The aim is to use the alignments to form a consensus sequence using \code{\link{consensusReadSeq}}.
This allows us to obtain accurate sequences of the original molecule, overcoming the high error rate of ONT.

The banded T-Coffee algorithm from \pkg{RSeqAn} is used to perform efficient multiple sequence alignment of many reads.
Note that, for DNA sequences, we assume that all reads are on the same strand.
This should automatically be the case if \code{\link{adaptorAlign}} was used with two different adaptors on either end of the read.
However, if the same adaptors were used on both ends, the strand of the sequence is unknown and may require other information to discern.

T-Coffee is not aware of quality strings, so we mask low-quality bases prior to the multiple sequence alignment.
All bases with error probabilities above \code{max.error} are set to \code{"N"}.
This ensures that poor base calls do not compromise the alignment quality.
If \code{min.qual=NA}, no masking is performed.
All masks are removed in the output strings unless \code{keep.mask=TRUE}, in which case the strings may contain \code{"N"}'s.
}

\value{
A \linkS4class{DataFrame} containing \code{alignments}, a \linkS4class{List} of character vectors.
Each vector corresponds to a read group in \code{groups} and contains the multiple sequence alignment strings for that group.

If \code{reads} is a QualityScaledDNAStringSet, the output DataFrame will also contain \code{qualities}, a List of quality score objects.
Each object contains the quality scores for each sequence in the MSA in the corresponding entry of \code{alignments}.
}

\author{
Florian Bieberich,
with modifications from Aaron Lun
}

\seealso{
    \code{\link{umiGroup}},
    \code{\link{consensusReadSeq}}
}

\references{
Reinert K, Dadi TH, Ehrhardt M \emph{et al.} (2017). 
The SeqAn C++ template library for efficient sequence analysis: A resource for programmers. 
\emph{J. Biotechnol.} 261:157-168

Notredame C, Higgins DG, Heringa J (2000). 
T-Coffee: A novel method for fast and accurate multiple sequence alignment. 
\emph{J. Mol. Biol.} 302(1):205-17
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


