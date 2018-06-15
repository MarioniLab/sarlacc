\name{expectedDist}
\alias{expectedDist}

\title{Calculate expected distances}
\description{Calculate expected distances between subsequences of the adaptor that should be identical across reads.}

\usage{
expectedDist(align.stats, position=NULL, offset=NULL, 
    number=100, get.seq=FALSE, max.err=NA)
}

\arguments{
\item{align.stats}{A DataFrame of alignment information, including the alignment strings for the read and adaptor.}
\item{position}{An integer vector of length 2 specifying the positions of first and last base of the UMI on the adaptor sequence.}
\item{offset}{An integer scalar in \code{[0, UMI length]} specifying the offset with which subsequences are extracted from the adaptor.}
\item{number}{An integer scalar specifying the number of reads to process.}
\item{get.seq}{A logical scalar indicating whether the extracted sequences should be returned.}
\item{max.err}{A numeric scalar specifying the maximum error probability above which bases will be masked.}
}

\details{
This function extracts a subsequence of the adaptor sequence from every read, and computes the distances between the extracted sequences from different reads. 
The aim is to provide an expectation for the distance for identical subsequences, given that all reads should originate from molecules with the same adaptor.
In this manner, we can obtain an appropriate threshold for \code{\link{umiGroup}} that accounts for sequencing and amplification errors.

By default, the extracted subsequence is taken from the intervals flanking the UMI region.
This ensures that the error rate in the extracted subsequence is as similar as possible to the UMI at that position on the read.
Specifically, for \code{position=c(START, END)} where \code{LEN=END-START+1}, 
we take the sequences at the intervals \code{c(START-offset, START-1)} and \code{c(END+1, END+LEN-offset)} for each read, and concatenate them together.

The function will infer \code{position} from the stretch of \code{N}'s in the adaptor sequence, though this can be explicitly specified as well.
By default, \code{offset} is set to half of \code{LEN}.
If \code{offset=0}, only the flanking interval on the 3' end of the UMI on the adaptor sequence is used.
Conversely, if \code{offset=LEN}, only the flanking interval on the 5' end is used.

Pairwise Levenshtein distances are computed between all extracted sequences.
This is quite computationally expensive, so we only process a random subset of these sequences by setting \code{number}.
If \code{align.stats} contains quality scores, bases with error probabilities above \code{max.qual} are replaced with \code{N}s.
Any \code{N}s are treated as missing and will contribute a mismatch score of 0.5, even for matches to other \code{N}s.
}

\value{
If \code{get.seq=FALSE}, a vector of pairwise distances between sequences that \emph{should} be identical.

Otherwise, a DNAStringSet or QualityScaledDNAStringSet object is returned containing the extracted sequence for each read.
}

\author{
Florian Bieberich,
with modifications by Aaron Lun
}

\seealso{
\code{\link{umiExtract}},
\code{\link{umiGroup}}
}

\examples{
adaptor <- "ACTACGACNNNNCAACGACT"
read <- c("ACTACGAGACTGCAACGACTCACAGCATCA",
          "ACTACGACACTGCAACGACTCACAGCATCA",
          "ACTACGATACTGCAACGACTCACAGCATCA",
          "ACTACGACACTGCAACGACTCACAGCATCA")

pw.aln <- pairwiseAlignment(pattern=read, subject=adaptor,
    type="local-global")
stats <- DataFrame(read=as.character(pattern(pw.aln)),
    adaptor=as.character(subject(pw.aln))) 

expectedDist(stats)
}
