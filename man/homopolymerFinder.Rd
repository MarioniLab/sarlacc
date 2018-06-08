\name{homopolymerFinder}
\alias{homopolymerFinder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find homopolymers}
\description{Identify homopolymers from the sequence. 
The definition of homopolymer is repeat of a single nucleotide consencutively.}
\usage{
homopolymerFinder(seq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq}{A \code{DNAStringSet} containing query sequence(s).}
}
\details{This function identifies regions of homopolymer in an input sequence and ignores gaps in the sequences. }
\value{An \code{IRangesList} object contianing regions, lengths and the base of homopolymers}
\author{}

\seealso{
\code{\link{homopolymerMatcher}}
}
\examples{
example_seq <- DNAStringSet("ATGGGGGTCGATATTTAAAAT")
homopolymerFinder(example_seq)
}