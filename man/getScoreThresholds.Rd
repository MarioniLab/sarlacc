\name{getScoreThresholds}
\alias{getScoreThresholds}

\title{Get score thresholds for adaptor alignments}
\description{Compute score thresholds to be used to determine if adaptors are present at the read ends.}

\usage{
getScoreThresholds(aligned, error=0.01, BPPARAM=SerialParam())
}

\arguments{
\item{aligned}{A list of adaptor alignment information and read sequences, identical to the output of \code{\link{adaptorAlign}}.}
\item{error}{A numeric scalar specifying the acceptable false discovery rate.}
\item{BPPARAM}{A BiocParallelParam object specifying whether alignment should be parallelized.
Currently only effective up to a maximum of 4 workers.}
} 

\value{
A named list with the following components:
\describe{
\item{\code{threshold1}:}{A numeric scalar specifying the alignment score threshold for adaptor 1.}
\item{\code{threshold2}:}{A numeric scalar specifying the alignment score threshold for adaptor 2.}
\item{\code{scores1}:}{A list containing \code{reads}, a sorted vector of alignment scores from \code{aligned$adaptor1};
and \code{scrambled}, a sorted vector of \dQuote{null} alignment scores for adaptor 1 against a scrambled read sequence.}
\item{\code{scores2}:}{Same as \code{scores1} but for adaptor 2.}
}
}

\details{
Read sequences (or specifically, the start and end of them, according to the tolerance specified in \code{\link{adaptorAlign}}) are scrambled,
and adaptor alignment is performed against the scrambled sequences.
For each adaptor, the scrambled alignment scores represent a null distribution of alignment scores.

We aim to define a score threshold above which the adaptor is considered to be present in the read sequence.
This is done by estimating the false discovery rate as the ratio of the number of scrambled scores above the threshold to the number of real alignment scores above the threshold.
We choose a threshold such that the false discovery rate is kept below the specified \code{error}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{adaptorAlign}},
\code{\link{chopReads}}
}

\examples{
example(adaptorAlign)
getScoreThresholds(out)
}
