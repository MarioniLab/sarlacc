\name{realizeReads}
\alias{realizeReads}

\title{Realize read sequences}
\description{Execute delayed operations (reverse-complementing, trimming) to generate the actual read sequences.}

\usage{
realizeReads(aligned, block.size=1e8, trim=TRUE) 
}

\arguments{
\item{aligned}{A \linkS4class{DataFrame} of alignment results produced by \code{\link{adaptorAlign}} - or, if \code{trim=TRUE}, by \code{\link{filterReads}}.}
\item{block.size}{Integer scalar specifying the number of bytes to read in at once from the FASTQ file, see \code{?\link{FastqStreamer}}.}
\item{trim}{Logical scalar indicating whether adaptor trimming should be performed.}
}

\details{
This function will parse the FASTQ file to extract the read sequences corresponding to the row names of \code{aligned}.
It will then perform reverse-complementing according to \code{aligned$reversed} to coerce all reads into the canonical orientation, see \code{?\link{adaptorAlign}} for details.
Finally, if \code{trim=TRUE}, it will trim off the adaptor sequences on the ends of the reads.

The reverse complementing and trimming can be envisaged as delayed operations from the previous functions (\code{\link{adaptorAlign}} and \code{\link{filterReads}}, respectively).
This sacrifices speed for memory efficiency in those steps, by avoiding the need to explicitly load the read sequences into memory.
Here, we have assumed that the number of reads is sufficiently reduced (e.g., by filtering and demultiplexing) to be able to load all reads.
}

\value{
A \linkS4class{QualityScaledDNAStringSet} object containing read sequences.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{adaptorAlign}} and \code{\link{filterReads}} for generation of \code{aligned}.
}

\examples{
example(filterReads)
realizeReads(filt)
}
