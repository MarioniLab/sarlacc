\name{groupReads}
\alias{groupReads}

\title{Group reads by sequence similarity}
\description{Group reads based on sequence similarity, using fast pairwise alignment from minimap2.}

\usage{
groupReads(reads, truncate=1000, tmp.dir=NULL, mm.cmd="minimap2", 
    mm.args=NULL, min.identity=0.7)
}

\arguments{
\item{reads}{A \code{DNAStringSet} object contain read sequences.}
\item{truncate}{An integer scalar specifying the length of sequence at the start of the read to use for alignment.
This can be set to \code{Inf} for no truncation, at the cost of more computational work.}
\item{working.dir}{A string containing the path to the directory to save intermediate files.
Defaults to a temporary directory in the working directory that is removed upon function exit.}
\item{mm.cmd}{A string containing the location to the minimap2 binary.}
\item{mm.args}{A character vector of further arguments to pass to minimap2.} 
\item{min.identity}{A numeric scalar specifying the minimum sequence identity between reads required for them to be assigned to the same cluster.}
}

\details{
Consider a graph where each node is a read.
Pairs of nodes are connected if they have sequence identity (PID4, see \code{?\link{pid}}) greater than \code{min.identity}.
Each component (i.e., one set of connected nodes) is defined as a group of reads.
For speed, we do not perform an exhaustive pairwise alignment but instead use the fast all-versus-all alignment mode in minimap2. 
}

\value{
A list of integer vectors.
Each vector corresponds to a group and contains the indices of \code{reads} in that single group.
}

\references{
Li, H. (2018).
Minimap2: pairwise alignment for nucleotide sequences. 
\emph{Bioinformatics}.
}

\seealso{
See \code{\link{umiGroup}} for the equivalent operation based on UMIs.

See \code{\link{multiReadAlign}} for the function that ultimately uses the grouping information.
}

\author{
Cheuk-Ting Law, 
with modifications by Aaron Lun.
}
