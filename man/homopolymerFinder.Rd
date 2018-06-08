\name{homopolymerFinder}
\alias{homopolymerFinder}

\title{Find homopolymers}
\description{Find homopolymer runs in a set of sequences.}

\usage{
homopolymerFinder(seq)
}

\arguments{
    \item{seq}{A DNAStringSet object.}
}

\details{
This function will identify homopolymers in a given set of sequences.
It is useful for investigating the homopolymer frequency in unknown sequences such as UMIs.
If the sequence is known, it is often more informative to use \code{\link{homopolymerMatcher}} instead.

Gapped sequences are supported - gaps will be ignored when considering homopolymer runs and computing coordinates.
However, ambiguous bases are not be given special treatment, and will be handled like any other IUPAC character.
}

\value{
An IRangesList object where each entry corresponds to a sequence in \code{seq}.
Each IRanges specifies the coordinates of a homopolymer run in the current sequence, along with the base for which it is a homopolymer.
}

\author{
Aaron Lun, with contributions from Cheuk-Ting Law
}

\seealso{
\code{\link{homopolymerMatcher}}
}

\examples{
seq <- DNAStringSet(c("AAAAAGGGGGCCCCCCTTTTT", "AAAAGGGGGCCCCTTTTT"))
homopolymerMatcher(seq)
}
