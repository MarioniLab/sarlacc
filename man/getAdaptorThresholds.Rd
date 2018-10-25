\name{getAdaptorThresholds}
\alias{getAdaptorThresholds}

\title{Get score thresholds for adaptor alignments}
\description{Compute score thresholds to be used to determine if adaptors are present at the read ends.}

\usage{
getAdaptorThresholds(aligned, error=0.01, number=1e5, BPPARAM=SerialParam())
}

\arguments{
\item{aligned}{A \linkS4class{DataFrame} of adaptor alignment information and read sequences, identical to the output of \code{\link{adaptorAlign}}.}
\item{error}{A numeric scalar specifying the acceptable false discovery rate.}
\item{number}{Integer scalar specifying the number of records to read at once from the FASTQ file, see \code{?\link{FastqStreamer}}.}
\item{BPPARAM}{A BiocParallelParam object specifying whether alignment should be parallelized.}
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

For the sake of speed, the function will only use a randomly sampled subset of reads (specified by \code{number}).
It is unnecessary to scramble every read in order to compute a satisfactorily precise threshold.
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
getAdaptorThresholds(out)
}
