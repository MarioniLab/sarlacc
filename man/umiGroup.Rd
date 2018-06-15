\name{umiGroup}
\alias{umiGroup}

\title{Group equivalent UMIs}
\description{Identify reads originating from the same DNA molecule, based on similar sequences for the unique molecular identifiers (UMIs).}

\usage{
umiGroup(UMI1, max.lev1=3, UMI2=NULL, max.lev2=max.lev1, max.err=NA)
}

\arguments{
\item{UMI1}{A DNAStringSet or QualityScaledDNAStringSet object containing the UMI sequence for each read.}
\item{max.lev1}{A integer scalar specifying the maximum distance between the UMI sequences of two reads in order for them to be considered equivalent.}
\item{UMI2}{An optional DNAStringSet or QualityScaledDNAStringSet object containing the second UMI sequence for each read.}
\item{max.lev2}{A integer scalar specifying the maximum distance between the UMI sequences of two reads in order for them to be considered equivalent.}
\item{max.err}{A numeric scalar specifying the maximum error probability, above which bases are to be masked.}
}

\details{
This function will identify groups of reads originating from the same DNA molecule, based on the presence of similar UMI sequences.
Each read is represented by a UMI sequence in \code{UMI1}, obtained using functions like \code{\link{umiExtract}}.

Similarity between UMI sequences is computed by \code{\link{stringDist}} using the Levenshtein distance.
Bases with error probabilities above \code{max.err} will be masked and replaced with \code{N}s, which contribute 0.5 mismatches to the distance calculation.
We then define UMI groups using a variation of the graph-based approach described by Smith et al. (2017):
\enumerate{
    \item Each UMI is a node, which is connected to all other UMIs with Levenshtein distances less than or equal to \code{max.lev1}.
    The number of edges is used as an estimate of the node density.
    \item Starting from the most dense node (denoted \code{source1}), we identify all nodes that are connected to it via paths of non-increasing density.
    In other words, the density along the path from the \code{source1} to the target node should not increase.
    These nodes are putatively assigned to \code{source1}.
    \item We repeat this with the next source node, i.e., the next most-dense node that has not already been assigned (denoted \code{source2}).
    Previously assigned nodes may be re-assigned to \code{source2} if it is connected via a shorter path than to \code{source1}.
    or, if the path lengths are equal, the immediate precursor node on the path back to \code{source2} is denser than the immediate precursor on the path back to \code{source1}.
    \item The above step is repeated until all nodes are assigned to a source node.
    The assignments define the groupings.
}
Grouping by non-increasing densities reflects the expected result from PCR decay products.

If \code{UMI2} is supplied, we assume that each read is represented by two UMI sequences (e.g., on both ends of the read).
Sequences in \code{UMI1} are compared to other sequences in \code{UMI1}, and sequences in \code{UMI2} are compared to other sequences in \code{UMI2}.
Reads are only connected in the graph if the distance for \code{UMI1} is less than \code{max.lev1} and the distance for \code{UMI2} is less than \code{max.lev2}.
}

\value{
An integer vector containing the groupings for all reads, where reads in the same group are assigned to the same level.
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
umiGroup(umi1, max.lev1=0)
   
# Now with qualities.
umiX <- QualityScaledDNAStringSet(umi1, 
    PhredQuality(c("11223344", "1123344", "12223344", "112233444")))
umiGroup(umiX, max.err=0.02)
umiGroup(umiX, max.err=0.02, max.lev1=1)

# Now with a UMI2.
umi2 <- DNAStringSet(c("AACCGGTT", "CACCGGTT", 
    "AACCCGGTA", "AACGGGTT"))
umiGroup(umi1, UMI2=umi2)
umiGroup(umi1, UMI2=umi2, max.lev2=0)
}
