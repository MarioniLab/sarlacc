\name{chopReads}
\alias{chopReads}

\title{Chop adaptors off reads}
\description{Trim the read sequences to remove adaptors where possible.}

\usage{
chopReads(aligned, score1, score2, essential1=TRUE, essential2=TRUE)
}

\arguments{
\item{aligned}{A list of adaptor alignment information and read sequences, identical to the output of \code{\link{adaptorAlign}}.}
\item{score1, score2}{Numeric scalars specifying the minimum alignment score required for positive match of adaptor 1 or 2 to the read sequence.}
\item{essential1, essential2}{Logical scalars indicating whether the presence of adaptor 1 or 2 is essential, i.e., the read should be discarded if this adaptor is not found.}
}

\details{
This function will chop off adaptor sequences so that they do not affect downstream analyses of the read sequences.
Identification of a matching adaptor sequence is based on the alignment scores exceeding a certain threshold, i.e., \code{score1} and \code{score2} for adaptors 1 and 2 respectively.
An appropriate threshold for each adaptor can be chosen by inspecting a histogram of scores across randomly chosen reads in the library.
This assumes that there is a clear distinction between the alignment qualities of reads with and without adaptors.
Ideally, two modes in the score histogram for a given adaptor would be observed, allowing users to set the threshold between the modes.

If \code{essential1=TRUE}, a read will be discarded if a positive match to adaptor 1 cannot be found.
This is useful in situations where the adaptor contains critical information such as the unique molecular identifier, which is required for downstream steps such as \code{\link{umiExtract}}.
However, if adaptor 1 is not important (e.g., only necessary for PCR amplification), then it does not matter that it cannot be identified.
In such cases, setting \code{essential1=FALSE} will ensure that reads without a match to adaptor 1 are not discarded (though the adaptor will still be trimmed if found).
The same logic applies to \code{essential2} for adaptor 2.
}

\value{
A list with the following components:
\describe{
\item{\code{reads}:}{A DNAStringSet containing the chopped read sequences. Reads without essential adaptors are discarded.}
\item{\code{adaptor1}:}{A DataFrame containing alignment strings for adaptor 1 to the original read sequence, only returned if \code{essential1=TRUE}.
Each entry corresponds to an element of the output \code{reads}.
This is similar to \code{aligned$adaptor1} except that all positional information on the pattern (i.e., read) has been removed, as this is no longer valid after trimming.}
\item{\code{adaptor2}:}{A DataFrame containing alignment strings for adaptor 2 to the original read sequence, only returned if \code{essential2=TRUE}.
Each entry corresponds to an element of \code{reads}.
This is similar to \code{aligned$adaptor2} except without positional information on the pattern.}
}
Other fields added by \code{\link{adaptorAlign}} are preserved and subsetted appropriately to correspond to the retained \code{reads}.
}

\author{
Florian Bieberich,
Aaron Lun
}

\seealso{
\code{\link{adaptorAlign}},
\code{\link{umiExtract}}
}

\examples{
example(adaptorAlign)

chopReads(out, score1=15, score2=20) # both identified

chopReads(out, score1=20, score2=20) # adaptor 1 not identified

chopReads(out, score1=20, score2=20, essential1=FALSE) # ... but that's okay.
}
