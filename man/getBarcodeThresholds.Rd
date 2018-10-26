\name{getBarcodeThresholds}
\alias{getBarcodeThresholds}

\title{Get score thresholds for barcode alignments}
\description{Compute score thresholds to be used to determine if the read sequence matches a known reference barcode.}

\usage{
getBarcodeThresholds(baligned, nmads=3)
}

\arguments{
\item{baligned}{A \linkS4class{DataFrame} of barcode alignment information for each read, identical to the output of \code{\link{barcodeAlign}}.}
\item{nmads}{A numeric scalar specifying the number of median absolute deviations to define outliers.}
} 

\value{
A named numeric vector of length 2.
This contains \code{score}, the threshold on the alignment score;
and \code{gap}, the threshold on the difference between the best and next-best score.
}

\details{
The threshold for each metric is defined at \code{nmads} MADs below the median across all reads.
The assumption here is that most reads are aligned to correct barcodes.
Any small outliers are considered to represent misalignments and are removed.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{barcodeAlign}} to generate \code{baligned}.
}

\examples{
example(barcodeAlign)
getBarcodeThresholds(debarcoded)
}
