\name{barcodeAlign}
\alias{barcodeAlign}

\title{Determine barcode identity}
\description{Determine the barcode identity for each read, based on alignment to a set of reference sequences.}

\usage{
barcodeAlign(align.stats, barcodes, position)
}

\arguments{
\item{align.stats}{A \linkS4class{DataFrame} of adaptor alignment information for all reads, equivalent to one of the \code{adaptor*} output fields from \code{\link{adaptorAlign}}.}
\item{barcodes}{A character vector or \linkS4class{DNAStringSet} containing the barcode sequences.}
\item{position}{An integer vector of length two specifying the start and end position of the barcode on the adaptor, equivalent to the same argument in \code{\link{umiExtract}}.}
}

\value{
A \linkS4class{DataFrame} containing:
\describe{
\item{\code{sequence}:}{A DNAStringSet object containing the read sequence corresponding to the putative barcode.}
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
\code{\link{adaptorAlign}} to generate \code{adaptor*} input. 

\code{\link{getBarcodeThresholds}} to determine a suitable threshold for \code{score}.
}

\examples{
example(adaptorAlign)

# Assigning reads to individual barcodes.
debarcoded <- barcodeAlign(out$adaptor1, 
    barcodes=c("AAAAAAAA", "TTTTTTTT", "CCCCCCCCC", "GGGGGGGG"),
    position=c(10, 17))
debarcoded
}
