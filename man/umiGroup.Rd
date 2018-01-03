\name{umiGroup}
\alias{umiGroup}

\title{Group equivalent UMIs}
\description{Identify reads originating from the same DNA molecule, based on similar sequences for the unique molecular identifiers (UMIs).}

\usage{
umiGroup(UMI1, threshold, UMI2=NULL, groups=NULL, flip=NULL) 
}

\arguments{
\item{UMI1}{A DNAStringSet or QualityScaledDNAStringSet object containing the UMI sequence for each read.}
\item{threshold}{A numeric scalar specifying the maximum distance between the UMI sequences of two reads in order for them to be considered equivalent.}
\item{UMI2}{A DNAStringSet or QualityScaledDNAStringSet object containing a second UMI sequence for each read, if present.}
\item{groups}{A factor specifying the groupings of reads, where only UMIs from the same group are compared to each other.}
\item{flip}{A logical vector indicating whether UMIs 1 and 2 should be flipped for each read.}
}

\details{
This function will identify groups of reads originating from the same DNA molecule, based on the presence of similar UMI sequences.
Each read is represented by a UMI sequence in \code{UMI1}, obtained using functions like \code{\link{umiExtract}}.

Similarity between UMI sequences is computed by \code{\link{stringDist}}, using the Levenshtein distance if \code{UMI1} is a DNAStringSet object.
Otherwise, the distance calculation will incorporate the quality information in the QualityScaledDNAStringSet object. 
Distances below \code{threshold} are used to define edges on an undirected graph where each read is represented by a node (Smith et al., 2017).
Each connected component of the graph defines a group of reads originating from the same molecule.

If \code{UMI2} is supplied, we assume that each read is represented by two UMI sequences (e.g., on both ends of the read).
The distance between a pair of reads is then defined as the sum of the distances for both UMIs.
Similarly, \code{threshold} applies to the sum of distances rather than the distance for each individual UMI.

By default, sequences in \code{UMI1} are compared to other sequences in \code{UMI1}, and sequences in \code{UMI2} are compared to other sequences in \code{UMI2}.
If \code{flip} is \code{TRUE} for any read, UMIs 1 and 2 for that read will be swapped before any comparison is performed.
This can be useful in some situations, e.g., to ensure that all read sequences are on the same strand if the adaptors are the same on both ends.

If \code{groups} is supplied, pairwise distances are only computed between pairs of reads in the same level of \code{groups}.
This is useful for reducing the computational complexity, as well as serving as a sanity check to only compare UMIs for reads that are actually similar.
Groupings can be defined in any manner; a simple approach is to use \code{\link{clusterReads}} and group reads based on similar alignment locations.

% In theory, we could identify which reads are on the opposite strand based on the UMIs alone.
% However, this results in some unresolvable edge cases whereby a group might consist of both the read *and* its reverse complement.
% This could arise if UMI1-UMI2 matches (directly or indirectly) to UMI2-UMI1, and it's not clear which one is right or wrong.
% It is safer to use the entire read sequence during alignment and grouping to determine the sequence, if that is necessary.
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
umiGroup(umi1, group=c(1,2,1,2))
umiGroup(umi1, threshold=0)
    
umi2 <- DNAStringSet(c("AACCGGTT", "CACCGGTT", 
    "AACCCGGTA", "AACGGGTT"))
umiGroup(umi1, umi2)
umiGroup(umi1, umi2, flip=c(TRUE, FALSE, FALSE, TRUE))
}
