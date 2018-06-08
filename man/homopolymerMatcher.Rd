\name{homopolymerMatcher}
\alias{homopolymerMatcher}

\title{Match homopolymer lengths}
\description{Match the true and observed homopolymer lengths based on pairwise alignments of read sequences to a reference.}

\usage{
homopolymerMatcher(alignments) 
}

\arguments{
    \item{alignments}{A Global PairwiseAlignmentsSingleSubject object, usually produced by \code{\link{pairwiseAlignment}}.
    The subject should be a constant reference sequence.}
}

\details{
This function will identify all \dQuote{true} homopolymers in the reference sequence (i.e., runs of the same base).
For each true homopolymer, it will identify the corresponding read subsequence based on the pairwise alignment.
The observed length of the true homopolymer is defined as the longest contiguous run of the same base (ignoring deletion characters) in the corresponding read subsequence.

This is most easily illustrated with a few examples below.
For demonstration purposes, only the true homopolymer region and the corresponding read subsequence are shown in uppercase.

\preformatted{  # Observed length of 2.
  acgtAA--tgca # Read
  acgtAAAAtgca # Reference

  # Observed length of 3, despite the deletion character.
  acgtAA-Atgca # Read
  acgtAAAAtgca # Reference

  # Observed length of 2, as the T breaks the run of A's.
  acgtAATAtgca # Read
  acgtAAAAtgca # Reference

  # Observed length of 3, before the breaking T.
  acgtAAATAtgca # Read
  acgt-AAAAtgca # Reference

  # Observed length of 6, including the insertions before and after.
  acgtAAAAAAtgca # Read
  acgt-AAAA-tgca # Reference

  # Observed length of 4, as the observed run must overlap actual homopolymer bases.
  acgtAAAATAAAAAAAAAtgca # Read
  acgtAAAA----------tgca # Reference
}
}

\value{
An IRanges object where each entry represents a homopolymer run in the reference sequence.
The metadata will contain \code{base}, the base identity of the run;
and \code{observed}, an IntegerList containing a named integer vector for each run.
Each integer vector contains the observed homopolymer lengths in the read sequences as the names,
and the frequency of those lengths as the integer values.
}

\author{
Aaron Lun, with contributions from Cheuk-Ting Law
}

\seealso{
\code{\link{pairwiseAlignment}},
\code{\link{homopolymerFinder}}
}

\examples{
aln <- pairwiseAlignment(subject=DNAString(c("AAAAGGGGGCCCCTTTT")), 
    DNAStringSet(c("AAAAAGGGGGCCCCCCTTTTT", "AAAAGGGGGCCCCTTTTT")))
homopolymerMatcher(aln)
}
