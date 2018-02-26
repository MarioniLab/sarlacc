\name{clusterReads}
\alias{clusterReads}

\title{Cluster reads}
\description{Cluster reads together based on whether they share sufficient identity after a pairwise alignment.}

\usage{
clusterReads(reads, groups, min.identity=80, 
    gapOpening=5, gapExtension=2, match=1, mismatch=0)
}

\arguments{
\item{reads}{A DNAStringSet or QualityScaledDNAStringSet object.}
\item{groups}{A factor specifying the groups in which each entry of \code{reads} belongs.}
\item{min.identity}{A numeric scalar specifying the minimum percentage of sequence identity, computed using \code{\link{pid}} with \code{type="PID3"}.}
\item{gapOpening, gapExtension}{Alignment parameters to be passed to \code{\link{pairwiseAlignment}}.}
\item{match, mismatch}{Alignment parameters to be passed to \code{\link{nucleotideSubstitutionMatrix}}, ignored if \code{reads} is a QualityScaledDNAStringSet object.}
}

\details{
For each group of reads, this function performs pairwise alignments and identifies clusters of reads with percentage identity no less than \code{min.identity} amongst its members.
Clustering is done in a single-linkage fashion where a read is included in a cluster if it has sufficient sequence identity to any existing read in the cluster.

This function relies on the fact that reads have already been grouped to some extent using UMIs in \code{\link{umiGroup}}. 
Further clustering within each group ensures that reads with very different sequences but similar UMIs (by chance or due to saturation of the UMI space) do not get aligned together.
}

\value{
An integer vector specifying the new cluster of each read.
This is guaranteed to be nested within \code{groups}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{umiGroup}},
\code{\link{tuneAlignment}}
}

\examples{
reads <- DNAStringSet(c(
    "AACGGGTCGCTCGACGAACGTCCGATGCAGCACG",
    "ACGGGTCGCTCGACGACTTCGTCGGACACTGCAGCACG",
    "AACGGGTCGCTCGAACTGCGTCCGATGCAGCACG",
    "AACGGGTCGCTCGATACGTCCGATGCAGCACG",
    "AACGGGTCGGTCGACGACTACGTCCGATGCAGCACG",
    "CGTGCGGTCATCGGACGTAGTCGTCGAGCGACCCGTT",
    "CGTGAGTCATCGGACGTAGTCGTCGACACCCCGTT",
    "CGTGCTAATCGGACGTAGTCGTCGAGCGCCCGTT"))
clusterReads(reads, integer(length(reads)))
}
