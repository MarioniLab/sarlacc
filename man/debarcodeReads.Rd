\name{debarcodeReads}
\alias{debarcodeReads}

\title{Determine barcode identity}
\description{Determine the barcode identity for each read, based on alignment to a set of reference sequences.}

\usage{
debarcodeReads(align.stats, barcodes, position)
}

\arguments{
\item{aligned}{A list of adaptor alignment information and read sequences, identical to the \code{adaptor*} output of \code{\link{adaptorAlign}}.}
\item{barcodes}{A character vector or \linkS4class{DNAStringSet} containing the barcode sequences.}
\item{position}{An integer vector of length two specifying the start and end position of the barcode on the adaptor, equivalent to the same argument in \code{\link{umiExtract}}.}
}

\value{
A \linkS4class{DataFrame} containing \code{id}, the index of the barcode in \code{barcodes} to which each read is assigned;
and \code{score}, the alignment score to the assigned barcode.
Each row corresponds to an entry in \code{align.stats}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{adaptorAlign}}
}

\examples{
example(adaptorAlign)
debarcodeReads(out$adaptor1,
    c("CGGG", "AAAT"),
    position=c(3, 7))
}
