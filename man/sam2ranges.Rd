\name{sam2ranges}
\alias{sam2ranges}

\title{SAM to genomic ranges}
\description{Import alignment positions in a SAM file to a GenomicRanges object.}

\usage{
sam2ranges(sam, minq = 10, restricted = NULL)
}

\arguments{
\item{sam}{String containing the path to a SAM file.}
\item{minq}{A numeric scalar specifying the minimum mapping quality score for an alignment to be retained.}
\item{restricted}{A character vector containing reference sequences to restrict to.}
}

\details{
We use SAM files as ONT CIGAR strings are too long to be handled with the BAM format, and \pkg{Rsamtools} is yet to support CRAM files.

If \code{restricted} is not \code{NULL}, only alignments to reference sequences in \code{restricted} will be considered.
This is useful for discarding alignments to scaffolds or other uninteresting parts of the genome.

If \code{minq=NULL}, no filtering is performed.
}

\value{
A GRangesList is returned where each entry corresponds to a read.
An entry may be of length greater than 1 if a read aligns to multiple locations (due to multi-mapping or chimerism).
}

\seealso{
\code{\link{GRangesList-class}},
\code{\link{cigarWidthAlongReferenceSpace}}
}

\author{
Florian Bieberich, with modifications from Aaron Lun
}

\examples{
library(Rsamtools)
fl <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
out <- asSam(fl, tempfile(fileext=""))
sam2ranges(out)

sam2ranges(out, restricted="seq1")

sam2ranges(out, minq=80)
}
