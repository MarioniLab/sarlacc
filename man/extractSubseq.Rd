\name{extractSubseq}
\alias{extractSubseq}

\title{Extract subsequence}
\description{Extract an arbitrary read subsequence corresponding to positions of the aligned adaptor.}

\usage{
extractSubseq(aligned, subseq1, subseq2, number=1e5, BPPARAM=SerialParam()) 
}

\arguments{
\item{aligned}{A \linkS4class{DataFrame} containing the output of \code{\link{adaptorAlign}}.}
\item{subseq1}{A list of two integer vectors \code{start} and \code{end} of equal length.
Parallel entries specify the start and end positions on adaptor 1 to extract the aligned read subequence.}
\item{subseq2}{Same as \code{subseq1} but for adaptor 2.}
\item{number}{Integer scalar specifying the number of records to read at once from the FASTQ file, see \code{?\link{FastqStreamer}}.}
\item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying how the parallelization is to be performed.}
}

\details{
This function will align the adaptors in \code{aligned} to the start and end of each read (see \code{?\link{adaptorAlign}}).
From the alignment, it will extract the subsequence of the read corresponding to the specified positions on the adaptor sequence in \code{subseq1} or \code{subseq2}.
This is useful in other functions such as \code{\link{expectedDist}}, which rely on read sequences corresponding to constant regions of the adaptor.

At least one of \code{subseq1} or \code{subseq2} must be specified.
}

\value{
A list containing up to two DataFrames.
Each DataFrame corresponds to an adaptor and contains the extracted read subsequences where each row corresponds to a row of \code{aligned}.
DataFrames are only returned for adaptors where \code{subseq*} was specified.
}

\seealso{
\code{\link{adaptorAlign}} to generate \code{aligned}.
}

\author{
Aaron Lun
}

\examples{
example(adaptorAlign)

# Let's say we want to take the first part of 'a1'.
substr(a1, 1, 9)
extractSubseq(out, subseq1=list(starts=1, ends=9))

# Let's say we also want to take some part of 'a2'.
substr(a2, 5, 11)
extractSubseq(out, subseq1=list(starts=1, ends=9),
    subseq2=list(starts=5, ends=11))
}
