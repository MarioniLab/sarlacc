\name{readsMerger}
\alias{readsMerger}

\title{Create consensus sequences}
\description{Group reads by using the pairwise alignment from minimap2 and UMI sequences. The grouped reads are merged together to create a consensus sequence.}

\usage{
    readsMerger(reads, umi, minimap2_dir = "", minimap2_arg = "-t 6", min_match = 0.7)
}

\arguments{
    \item{reads}{A \code{QualityScaledBStringSet} object containing reads and sequences quality.
    It is usually in fastq format and can be imported by \code{readQualityScaledDNAStringSet}}
    \item{umi}{A \code{DNAStringSet} object containing UMI sequnece and having the same order as the reads.}
    \item{minimap2_dir}{A directory to minimap2}
    \item{minimap2_arg}{Arugments used for minimap2 pairwise alignment. Please see (\link{https://github.com/lh3/minimap2})}
    \item{min_match}{Similarity (PID4) threshold for clustering reads.}
}
\details{
Pairs of similar reads are found by minimap2 pairwise alignment. 
The read pairs are only retained if they have a similarity (PID4) higher than \code{min_match}. 
For those remaining reads, \code{igraph} retrieves all reads that are connected and cluster them into groups.
Some sequences are not able to find any sequence to pair up due to their low expression and errors caused by Nanopore sequencing.
These sequences are kept for later analysis because although these sequences deviate too far from its replicates, they are still possible to be paired up once better consensus sequences appear in the later iterations.

After clustering reads by similarity, \code{\link{umiGroup2}} further subgroup UMIs in each cluster according to their Levenshtein distances.
Those subgrouped reads are considered to be originated from the same molecule during PCR amplification.
\code{\link{multiReadAlign}} executes multiple sequence alignment to the subgrouped reads and passes the result to \code{\link{consensusReadSeq}} for consensus sequence calling. 
In this step, the original sequence IDs of each consensus sequence are recorded into the \code{elementMetadata}.
    
Repeatedly, minimap2 carries out pairwise alignment to the consensus sequences and \code{\link{umiGroup2}} subgroups the reads in each cluster. 
Then a consensus sequence is created by using the previously stored sequences and the newly grouped sequences. 
The process of pairwise alignment, UMI grouping and consensus sequence calling are repeated until no more subgroup is found.
}


\value{
A \code{QualityScaledDNAStringSet} object containing merged consenus sequence and the sequence quality.
The detail of computed quality can be found in \code{\link{consensusReadSeq}}. 
Furthermore, the history of each consensus sequence is stored in the \code{elementMetadata} and the merged sequences are named by combining all their original sequence IDs and separated by a space.
}

\references{
Li, H. (2018).
Minimap2: pairwise alignment for nucleotide sequences. 
\emph{Bioinformatics}.
}