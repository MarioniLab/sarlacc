\name{groupAligned}
\alias{groupAligned}

\title{Group aligned reads by location}
\description{Group alignments based on their location in the reference, with an additional threshold on the minimum overlap.}

\usage{
groupAligned(locations, threshold=0.8, cluster.fun=cluster_fast_greedy)
}

\arguments{
\item{locations}{A \linkS4class{GenomicRanges} object containing an alignment location for each read.
Alternatively, a \linkS4class{GenomicRangesList} object containing a multi-entry alignment for each read.}
\item{threshold}{A numeric scalar specifying the minimum percent overlap between alignment locations.}
\item{cluster.fun}{A clustering function that takes a \linkS4class{graph} object and returns a \linkS4class{communities} object.}
}

\details{
Consider a graph where each node is a read.
Pairs of nodes are connected if their alignment widths overlap by at least \code{threshold}.
Communities of highly connected reads are then obtained using \code{\link{cluster.fun}}, where each community represents a group of similar reads.

The overlap proportion for each pair is based on the width of the overlapping region divided by the larger of the two alignment widths.
This means that the overlap region must be large relative to both reads, consistent with their origin from the same underlying molecule.
A multi-entry alignment for each read is supported assuming the various entries for a given read are non-overlapping.
}

\value{
An integer vector specifying the group identity for each entry of \code{locations}.
Unmapped reads have identities of \code{NA}.
}

\seealso{
See \code{\link{umiGroup}} for the equivalent operation based on UMIs.
}

\author{
Aaron Lun.
}

\examples{
library(Rsamtools)
fl <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
out <- asSam(fl, tempfile(fileext=""))
gr <- sam2ranges(out)

groups <- groupAligned(gr)
table(groups) # number of reads in each group.
}
