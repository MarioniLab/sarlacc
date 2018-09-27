\name{getBarcodeThresholds}
\alias{getBarcodeThresholds}

\title{Get score thresholds for barcode alignments}
\description{Compute score thresholds to be used to determine if the read sequence matches a known reference barcode.}

\usage{
getBarcodeThresholds(baligned, error=0.01, BPPARAM=SerialParam())
}

\arguments{
\item{baligned}{A \linkS4class{DataFrame} of barcode alignment information for each read, identical to the output of \code{\link{debarcodeReads}}.}
\item{error}{A numeric scalar specifying the acceptable false discovery rate.}
\item{BPPARAM}{A BiocParallelParam object specifying whether alignment should be parallelized.}
} 

\value{
A named numeric vector of length 2.
This contains \code{score}, the threshold on the alignment score;
and \code{gap}, the threshold on the difference between the best and next-best score.
}

\details{
Read sequences at the barcode location are scrambled, and barcode alignment is performed against the scrambled sequences.
We define a threshold on the scores above which the barcode is considered to be present, based on controlling the false discovery rate below \code{error} 
(see \code{\link{getAdaptorThreshold}} for details).

We also define a threshold on the \emph{difference} between the best and next-best score for each read.
It may be preferable to filter reads on the difference threshold if some barcodes are similar to each other.
Otherwise, the origin of a read is not clear if it can be feasibly assigned to multiple barcodes.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{adaptorAlign}},
\code{\link{chopReads}}
}

\examples{
example(debarcodeReads)
getBarcodeThresholds(debarcoded)
}
