\name{umiGroup}
\alias{umiGroup}

\title{Group equivalent UMIs}
\description{Identify reads originating from the same DNA molecule, based on similar sequences for the unique molecular identifiers (UMIs).}

\usage{
umiGroup(UMI1, threshold1=3, UMI2=NULL, threshold2=threshold1, max.err=NA, 
    groups=NULL)
}

\arguments{
\item{UMI1}{A \linkS4class{DNAStringSet} or \linkS4class{QualityScaledDNAStringSet} object containing the UMI sequence for each read.}
\item{threshold1}{A integer scalar specifying the maximum distance between the UMI sequences of two reads in order for them to be considered equivalent.}
\item{UMI2}{An optional DNAStringSet or QualityScaledDNAStringSet object containing the second UMI sequence for each read.}
\item{threshold2}{A integer scalar specifying the maximum distance between the UMI sequences of two reads in order for them to be considered equivalent.}
\item{max.err}{A numeric scalar specifying the maximum error probability above which bases are to be masked.}
\item{groups}{A factor specifying the group to which each read belongs.}
}

\details{
This function will identify groups of reads originating from the same DNA molecule, based on the presence of similar UMI sequences.
Each read is represented by a UMI sequence in \code{UMI1}, obtained from the read sequence using \code{\link{adaptorAlign}}.

We use a greedy approach to group reads from the same molecule.
For each UMI, we define neighboring UMIs as those within a Levenshtein distance of \code{threshold1}.
The first group is defined as all neighbors of the UMI with the most neighbors.
The process is then repeated with the remaining UMIs.
This is a conservative procedure that errs on the side of under-merging to avoid grouping reads from distinct molecules.

Bases with error probabilities above \code{max.err} will be masked and replaced with \code{N}s, which contribute 0.5 mismatches to the distance calculation.
If \code{UMI2} is supplied, we assume that each read is represented by two UMI sequences (e.g., on both ends of the read).
Sequences in \code{UMI1} are compared to other sequences in \code{UMI1}, and sequences in \code{UMI2} are compared to other sequences in \code{UMI2}.
UMIs are only defined as neighbors if the Levenshtein distances for the first and second UMIs are below \code{threshold1} and \code{threshold2}, respectively.

If \code{groups} is specified, reads are only grouped together if they have the same level of \code{groups}.
This is useful for restricting the UMI comparisons to pre-defined groups of reads, e.g., reads that map to the same transcript family.
}

\value{
A list of integer vectors, where each vector represents a group and contains the indices of the reads (i.e., with respect to \code{UMI1}) in that group.
}

\seealso{
\code{\link{adaptorAlign}} to obtain the UMI sequences.
}

\examples{
umi1 <- DNAStringSet(c("AACCGGTT", "AACGGTT", 
    "ACCCGGTT", "AACCGGTTT"))

umiGroup(umi1)
umiGroup(umi1, threshold1=0)
   
# Now with qualities.
umiX <- QualityScaledDNAStringSet(umi1, 
    PhredQuality(c("11223344", "1123344", "12223344", "112233444")))
umiGroup(umiX, max.err=0.02)
umiGroup(umiX, max.err=0.02, threshold1=1)

# Now with a UMI2.
umi2 <- DNAStringSet(c("AACCGGTT", "CACCGGTT", 
    "AACCCGGTA", "AACGGGTT"))
umiGroup(umi1, UMI2=umi2)
umiGroup(umi1, UMI2=umi2, threshold2=0)
}
