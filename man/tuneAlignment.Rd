\name{tuneAlignment}
\alias{tuneAlignment}

\title{Tune alignment parameters}
\description{Tune parameters for adaptor alignment to maximize discriminative power compared to a randomized control.}

\usage{
tuneAlignment(adaptor1, adaptor2, reads, tolerance=100, 
    gapOp.range=c(4, 10), gapExt.range=c(1, 5), 
    match.range=c(1, 2), mismatch.range=c(-1, 0)) 
}

\arguments{
\item{adaptor1, adaptor2}{A string or DNAString object containing the 5'-to-3' sequences of the adaptors on each end of the read.}
\item{reads}{A DNAStringSet object containing the read sequences.
A QualityScaledDNAStringSet can also be supplied containing quality scores.}
\item{tolerance}{An integer scalar specifying the region at the ends of the reads to search for the adaptors.}
\item{gapOp.range}{An integer vector of length 2 specifying the boundaries of the grid search for the gap opening penalties.}
\item{gapExt.range}{An integer vector of length 2 specifying the boundaries of the grid search for the gap extension penalties.}
\item{match.range}{An integer vector of length 2 specifying the boundaries of the grid search for the match score.}
\item{mismatch.range}{An integer vector of length 2 specifying the boundaries of the grid search for the mismatch score.}
}

\details{
This function will align adaptors to the start and end of read sequences in the same manner as \code{\link{adaptorAlign}}.
It will then perform a grid search to identify the best parameters for alignment.
This is done by repeating the alignments for all possible combinations of integer gap opening or extension penalties.
If \code{reads} is a DNAStringSet, all combinations of the match and mismatch scores will also be explored.

To evaluate each parameter combination, we examine the distribution of combined alignment scores for all reads.
This represents the best adaptor alignment and is equivalent to the approach used in \code{\link{adaptorAlign}} to determine the read orientation.
The best parameter combiantion is which minimizes the overlap between the distribution of maximum alignment scores for \code{reads} and that of a scrambled control.
Obviously, we only look for combinations where the former distribution is shifted towards higher scores compared to the scrambled control.
}

\value{
A list containing \code{parameters}, itself a list containing the optimal values of all specified alignment parameters.
(If \code{reads} is a QualityScaledDNAStringSet, the match and mismatch parameters are \code{NA}).

The top-level list will also contain \code{scores}, another list containing numeric vectors of alignment scores for the reads and scrambled controls.
These scores are computed at the optimal values and are reported for diagnostic purposes.
}

\author{Aaron Lun
}

\seealso{
\code{\link{adaptorAlign}}
}

\examples{
tuneAlignment("CGTACGACGAT", "TCGAGCGTTAC", c(
    "CGTACGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGCTCGA",
    "CGTCGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGCTCGA",
    "CGTACGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGCCGA",
    "CGTACGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGCTCGA",
    "CGCTACGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGGTCGA"))
}
