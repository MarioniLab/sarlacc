\name{errorFinder}
\alias{errorFinder}

\title{Find sequencing errors}
\description{Find errors in one or more read sequences compared to a single reference sequence, based on pairwise alignments between them.}

\usage{
errorFinder(alignments) 
}

\arguments{
    \item{alignments}{A Global PairwiseAlignmentsSingleSubject object, usually produced by \code{\link{pairwiseAlignment}}.
    The subject should be a constant reference sequence.}
}

\details{
For each position in the reference sequence, this function records the observed frequency of each base across all aligned read sequences.
A transition matrix is computed describing the number of times a base in the reference sequence (rows) is observed as a base in the read sequence (column).

Insertions in the read sequence are assigned to the base position that they precede in the reference sequence.
So an alignment like:

\preformatted{  AAAATGGGG # Read
  AAA--GGGG # Reference
}

would yield an insertion of length 2 that is assigned to base position 4 on the reference sequence.
Insertions at the end of the sequence are assigned to an \dQuote{imaginary} base position at one plus the sequence length.

We also record the observed frequency of deletions for each position in the reference.
}

\value{
A list is returned containing a DataFrame in \code{full} and a transition matrix in \code{transition}.

In the \code{full} DataFrame, each row represents a position on the reference sequence.
The DataFrame contains the fields:
\describe{
    \item{\code{base}:}{Character, the base at the current position.}
    \item{\code{A},\code{C},\code{T},\code{G}:}{Integer, the frequency of observing A, C, T or G at the current position.}
    \item{\code{deletion}:}{Integer, the frequency of observing a deletion at the current position.}
    \item{\code{insertion}:}{IntegerList, where the entry for each position is a named integer vector. 
        This describes a distribution of insertions immediately preceding the current position.
        The name of each element is the insertion length while the integer element itself is the frequency of that length.
    }
}
The last row corresponds to a hypothetical one-past-the-end position, and is \code{NA} for all fields except for \code{insertion} (i.e., when an insertion occurs at the sequence end).
}

\author{
Aaron Lun, with contributions from Cheuk-Ting Law
}

\seealso{
\code{\link{pairwiseAlignment}}
}

\examples{
aln <- pairwiseAlignment(subject=DNAString(c("AAACGATCAGCTACGAACACT")), 
    DNAStringSet(c("AACGAGGGCCACCTAGGAAGAAT", "AACCAATCCAGCTACGCAACGACT")))
errorFinder(aln)
}
