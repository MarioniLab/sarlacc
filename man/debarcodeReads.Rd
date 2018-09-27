\name{debarcodeReads}
\alias{debarcodeReads}

\title{Determine barcode identity}
\description{Determine the barcode identity for each read, based on alignment to a set of reference sequences.}

\usage{
debarcodeReads(align.stats, barcodes, position)
}

\arguments{
\item{align.stats}{A \linkS4class{DataFrame} of adaptor alignment information and read sequences, identical to the \code{adaptor*} output of \code{\link{adaptorAlign}}.}
\item{barcodes}{A character vector or \linkS4class{DNAStringSet} containing the barcode sequences.}
\item{position}{An integer vector of length two specifying the start and end position of the barcode on the adaptor, equivalent to the same argument in \code{\link{umiExtract}}.}
}

\value{
A \linkS4class{DataFrame} containing:
\describe{
\item{\code{sequence}:}{A DNAStringSet object containing the read sequence corresponding to the putative barcode.}
\item{\code{barcode}:}{Integer, the index of the best-matching barcode in \code{barcodes} to which each read is assigned.}
\item{\code{score}:}{Numeric, the alignment score to the assigned barcode.}
\item{\code{alternative}:}{Numeric, the alignment score to the next best-matching barcode.}
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
a1 <- "AACGGGTCGNNNNNNNACGNNNNACGA" # first N stretch is a barcode.
a2 <- "CGTGCTGCATCG"

# Mocking up reads with different barcodes (same UMI for simplicity).
reads <- DNAStringSet(c(
    "AACGGGTCGAAAAAAAACGACGTACGACGACTGACTGACTAGCATTGCGATGCAGCACG",
    "AACGGGTCGAATTAATTCGACGTACGACGACTGACTGACTAGCATTGCGATGCAGCACG",
    "AACGGGTCGTTAATTAACGACGTACGACGACTGACTGACTAGCATTGCGATGCAGCACG",
    "AACGGGTCGAAGGAAGGCGACGTACGACGACTGACTGACTAGCATTGCGATGCAGCACG",
    "AACGGGTCGCCAACCAACGACGTACGACGACTGACTGACTAGCATTGCGATGCAGCACG",
    "AACGGGTCGCCGGCCGGCGACGTACGACGACTGACTGACTAGCATTGCGATGCAGCACG",
    "AACGGGTCGCCTTCCTTCGACGTACGACGACTGACTGACTAGCATTGCGATGCAGCACG"
))

(out <- adaptorAlign(adaptor1=a1, adaptor2=a2, reads=reads))

# Assigning reads to individual barcodes.
debarcoded <- debarcodeReads(out$adaptor1, 
    barcodes=c("AAAAAAAA", "AATTAATT", "TTAATTAA", "AAGGAAGG", 
        "CCAACCAA", "CCGGCCGG", "CCTTCCTT"),
     position=c(10, 17))
            
}
