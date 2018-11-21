\name{groupAligned}
\alias{groupAlignedGenome}
\alias{groupAlignedTrans}

\title{Group aligned reads by location}
\description{Group alignments based on their location in the reference, with an additional threshold on the minimum overlap.}

\usage{
groupAlignedGenome(locations)

groupAlignedTrans(locations)
    
}

\arguments{
\item{locations}{A \linkS4class{GenomicRanges} object containing an alignment location for each read.
Alternatively, a \linkS4class{GenomicRangesList} object containing a multi-entry alignment for each read.}
}

\details{
For \code{groupAlignedGenome}, aligned reads are grouped together based on overlaps to intervals on the reference genome.
We define each interval of interest as a contiguous region of unbroken coverage along a given strand of a reference sequence.
Each aligned read is then assigned to the first interval that it overlaps, and the set of reads assigned to each interval forms a group.

For \code{groupAlignedTrans}, aligned reads are grouped together according to the name and strand of the reference sequence to which it aligns.
This is faster and more relevant for alignments to transcriptomes containing relatively short reference sequences.
All reads assigned to a particular name/strand combination form a group.

For multi-entry alignments, each element of \code{locations} will be sorted prior to identifying groups.
This reduces the chance that reads with similar alignment locations are assigned to different groups solely due to the ordering of alignment entries.
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
