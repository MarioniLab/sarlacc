\name{umiExtract}
\alias{umiExtract}

\title{Extract UMIs}

\description{Extract unique molecular identifier (UMI) sequences from an alignment of a UMI-containing adaptor to the read sequence.}

\usage{
umiExtract(align.stats, position=NULL) 
}

\arguments{
\item{align.stats}{A DataFrame of alignment information, including the alignment strings for the read and adaptor.}
\item{position}{An integer vector of length 2 specifying the positions of first and last base of the UMI on the adaptor sequence.}
}

\details{
This function extracts the UMI sequence from the alignment between the adaptor and UMI sequences.
The alignment should be global with respect to the adaptor sequence, i.e., starting at position 1.
The location of the UMI within the adaptor is defined by \code{position} using 1-based coordinates. 
If \code{postion} is not supplied, the location of the UMI is automatically determined based on the location of \code{N}s in the adaptor sequence.
}

\value{
A DNAStringSet of UMI sequences, one for each row of \code{align.stats}.
}

\seealso{
\code{\link{adaptorAlign}}
}

\author{
Florian Bieberich,
Aaron Lun
}

\examples{
adaptor <- "ACTACGACNNNNCAACGACT"
read <- "ACTACGACACTGCAACGACTCACAGCATCA"

# Mocking up an adaptor alignment: use adaptorAlign() for real data.
pw.aln <- pairwiseAlignment(pattern=read, subject=adaptor,
    type="local-global")
stats <- DataFrame(read=as.character(pattern(pw.aln)),
    adaptor=as.character(subject(pw.aln))) 
rownames(stats) <- "WHEE"

(out <- umiExtract(stats))

# Adding some qualities.
stats$quality <- PhredQuality("983475298471982348712367432093")
(out <- umiExtract(stats))
}
