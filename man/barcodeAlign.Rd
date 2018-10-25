\name{barcodeAlign}
\alias{barcodeAlign}

\title{Determine barcode identity}
\description{Determine the barcode identity for each read, based on alignment to a set of reference sequences.}

\usage{
barcodeAlign(sequences, barcodes, gapOpening=5, gapExtension=1, 
    BPPARAM=SerialParam())
}

\arguments{
\item{sequences}{A \linkS4class{QualityScaledDNAStringSet} of observed barcode sequences, usually one of the subsequences from \code{\link{adaptorAlign}}.}
\item{barcodes}{A character vector or \linkS4class{DNAStringSet} containing the reference barcode sequences.}
\item{gapOpening}{A numeric scalar specifying the gap opening penalty.}
\item{gapExtension}{A numeric scalar specifying the gap extension penalty.}
\item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying whether alignment should be parallelized.}
}

\details{
Each entry of \code{sequences} is expected to contain the read subsequence that was aligned to the (usually ambiguous) barcode region of the adaptor in \code{\link{adaptorAlign}}.
These represent observed barcodes that are aligned to the reference set of known barcode sequences.
The read corresponding to each entry of \code{sequences} is assigned to the reference barcode to which the subsequence aligns best.
}

\value{
A \linkS4class{DataFrame} containing:
\describe{
\item{\code{barcode}:}{Integer, the index of the best-matching barcode in \code{barcodes} to which each read is assigned.}
\item{\code{score}:}{Numeric, the alignment score to the assigned barcode.}
\item{\code{gap}:}{Numeric, the difference of \code{score} from the alignment score to the next best-matching barcode.}
}
Each row corresponds to a read (i.e., row) in \code{align.stats}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{adaptorAlign}} to generate subsequences.

\code{\link{getBarcodeThresholds}} to determine a suitable threshold for \code{score}.
}

\examples{
example(adaptorAlign)

# Assigning reads to individual barcodes.
debarcoded <- barcodeAlign(out$adaptor1$subseqs$Sub1, 
    barcodes=c("AAAAAAAA", "TTTTTTTT", "CCCCCCCCC", "GGGGGGGG"),
    gapOpening=metadata(out$adaptor1)$gapOpening,
    gapExtension=metadata(out$adaptor1)$gapExtension)
debarcoded
}
