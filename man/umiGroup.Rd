\name{umiGroup}
\alias{umiGroup}

\title{Group equivalent UMIs}
\description{Identify reads originating from the same DNA molecule, based on similar sequences for the unique molecular identifiers (UMIs).}

\usage{
umiGroup(UMI1, UMI2=NULL, groups=NULL, threshold=2, flip=FALSE) 
}

\arguments{
\item{UMI1, UMI2}{A DNAStringSet object containing the UMI sequence(s) for each read.}
\item{groups}{A factor specifying the groupings of reads, where only UMIs from the same group are compared to each other.}
\item{threshold}{An integer scalar specifying the maximum Levenshtein distance between the UMI sequences of two reads in order for them to be considered equivalent.}
\item{flip}{A logical scalar indicating whether strandedness should be ignored.}
}

\details{
This function will identify groups of reads originating from the same DNA molecule, based on the presence of similar UMI sequences.
Each read is represented by a UMI sequence in \code{UMI1}, obtained using functions like \code{\link{umiExtract}}.

By default, the Levenshtein distance between all pairs of reads in \code{UMI1} is computed.
Distances below \code{threshold} are used to define edges on an undirected graph where each read is represented by a node (Smith et al., 2017).
Each connected component of the graph defines a group of reads originating from the same molecule.

If \code{UMI2} is supplied, we assume that each read is represented by two UMI sequences (e.g., on both ends of the read).
The Levenshtein distance between a pair of reads is then defined as the sum of the distances for both UMIs.
Similarly, \code{threshold} applies to the sum of distances rather than the distance for each individual UMI.

By default, sequences in \code{UMI1} are compared to other sequences in \code{UMI1}, and sequences in \code{UMI2} are compared to other sequences in \code{UMI2}.
If \code{flip=TRUE}, sequences in \code{UMI1} will be compared to sequences in \code{UMI2} and vice versa.
This can be useful in some situations, e.g., when each UMI is obtained from identical adaptors at both end of the reads.

If \code{groups} is supplied, pairwise distances are only computed between pairs of reads in the same level of \code{groups}.
This is useful for reducing the computational complexity, as well as serving as a sanity check to only compare UMIs for reads that are actually similar.
Groupings can be defined in any manner; a simple approach is to use \code{\link{clusterReads}} and group reads based on similar alignment locations.
}

\value{
An integer vector containing the groupings for all reads, where reads in the same group have the same level.
}

\seealso{
\code{\link{umiExtract}},
\code{\link{clusterReads}}
}

\references{
Smith T, Heger A, Sudbery I (2017). 
UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. 
\emph{Genome Res.}, 27, 3:491-499.
}

\examples{
umi1 <- DNAStringSet(c("AACCGGTT", "AACGGTT", 
    "ACCCGGTT", "AACCGGTTT"))

umiGroup(umi1)
umiGroup(umi1, group=c(1,2,1,2))
umiGroup(umi1, threshold=0)
    
umi2 <- DNAStringSet(c("AACCGGTT", "CACCGGTT", 
    "AACCCGGTA", "AACGGGTT"))
umiGroup(umi1, umi2)
umiGroup(umi1, umi2, flip=TRUE)
}
