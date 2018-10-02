\name{tuneAlignment}
\alias{tuneAlignment}

\title{Tune alignment parameters}
\description{Tune parameters for adaptor alignment to maximize discriminative power compared to a randomized control.}

\usage{
tuneAlignment(adaptor1, adaptor2, filepath, tolerance=200, 
    number=10000, gapOp.range=c(4, 10), gapExt.range=c(1, 5), 
    qual.type=c("phred", "solexa", "illumina"), BPPARAM=SerialParam()) 
}

\arguments{
\item{adaptor1, adaptor2}{A string or DNAString object containing the 5'-to-3' sequences of the adaptors on each end of the read.}
\item{filepath}{A string containing the path to the FASTQ file, or a connection object to a FASTQ file.}
\item{tolerance}{An integer scalar specifying the length of the ends of the reads to search for adaptors.}
\item{number}{An integer scalar specifying the number of randomly sampled reads to use for tuning.}
\item{gapOp.range}{An integer vector of length 2 specifying the boundaries of the grid search for the gap opening penalties.}
\item{gapExt.range}{An integer vector of length 2 specifying the boundaries of the grid search for the gap extension penalties.}
\item{qual.type}{String specifying the type of quality scores in \code{filepath}.}
\item{BPPARAM}{A BiocParallelParam object specifying whether alignment should be parallelized.
Currently only effective up to a maximum of 4 workers.}
}

\details{
This function will align adaptors to the start and end of read sequences in the same manner as \code{\link{adaptorAlign}}.
It will then perform a grid search to identify the best parameters for alignment.
This is done by repeating the alignments for all possible combinations of integer gap opening or extension penalties.

To evaluate each parameter combination, we examine the distribution of combined alignment scores for all reads.
This represents the best adaptor alignment and is equivalent to the approach used in \code{\link{adaptorAlign}} to determine the read orientation.
The best parameter combiantion is which minimizes the overlap between the distribution of maximum alignment scores for \code{reads} and that of a scrambled control.
Obviously, we only look for combinations where the former distribution is shifted towards higher scores compared to the scrambled control.
}

\value{
A list containing \code{parameters}, itself a list containing the optimal values of all specified alignment parameters.
The top-level list will also contain \code{scores}, another list containing numeric vectors of alignment scores for the reads and scrambled controls at the optimal parameters.
}

\author{Aaron Lun
}

\seealso{
\code{\link{adaptorAlign}} to use these parameters.
}

\examples{
# Mocking up a small data set.
a1 <- "AACGGGTCGNNNNNNNACGTACGTNNNNACGA" 
a2 <- "CGTGCTGCATCG"
fout <- tempfile(fileext=".fastq")
ref <- sarlacc:::mockReads(a1, a2, fout, nmolecules=1, 
    nreads.range=c(10, 10), seqlen.range=c(50, 200))

# Aligning it.
(out <- tuneAlignment(adaptor1=a1, adaptor2=a2, filepath=fout))
}
