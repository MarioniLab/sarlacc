\name{filterReads}
\alias{filterReads}

\title{Filter out reads without adaptors}
\description{Filter reads to only retain those with the essential adaptors.}

\usage{
filterReads(aligned, score1, score2, essential1=TRUE, essential2=TRUE)
}

\arguments{
\item{aligned}{A list of adaptor alignment information and read sequences, identical to the output of \code{\link{adaptorAlign}}.}
\item{score1, score2}{Numeric scalars specifying the minimum alignment score required for positive match of adaptor 1 or 2 to the read sequence.}
\item{essential1, essential2}{Logical scalars indicating whether the presence of adaptor 1 or 2 is essential, i.e., the read should be discarded if this adaptor is not found.}
}

\details{
Identification of a matching adaptor sequence is based on the alignment scores exceeding a certain threshold, i.e., \code{score1} and \code{score2} for adaptors 1 and 2 respectively.
An appropriate threshold for each adaptor can be chosen with methods like \code{\link{getAdaptorThresholds}}.

If \code{essential1=TRUE}, a read will be discarded if a positive match to adaptor 1 cannot be found.
This is useful when the adaptor contains critical information such as the unique molecular identifier, which is required for downstream steps such as \code{\link{umiExtract}}.
However, if adaptor 1 is not important (e.g., only necessary for PCR amplification), then it does not matter that it cannot be identified.
In such cases, setting \code{essential1=FALSE} will ensure that reads without a match to adaptor 1 are not discarded.
The same logic applies to \code{essential2} for adaptor 2.

This function will also report the start and end of the read sequence between the adaptors, marking the interval that would remain after adaptor trimming.
These positions are reported relative to the canonical orientation, see \code{?\link{adaptorAlign}} for more details.
The interval may include the start or end of the read if adaptor 1 or 2, respectively, is not essential and not found. 
Reads will also be discarded if adaptor 1 ends after adaptor 2 begins, as this implies that the adaptors overlap (and that there is no sequence in between).
}

\value{
\code{aligned} is returned, subsetted to only retain rows (i.e., reads) that have all essential adaptors.
It also contains the additional \code{trim.start} and \code{trim.end} fields, indicating the start of the sequence between the two adaptors.
}

\author{
Florian Bieberich,
Aaron Lun
}

\seealso{
\code{\link{adaptorAlign}} for the input into this function.
}

\examples{
example(adaptorAlign)

(filt <- filterReads(out, score1=5, score2=5)) # both identified

filterReads(out, score1=100, score2=5) # adaptor 1 not identified

filterReads(out, score1=100, score2=5, essential1=FALSE) # ... but that's okay.
}
