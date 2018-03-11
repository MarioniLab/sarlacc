\name{adaptorAlign}
\alias{adaptorAlign}

\title{Align adaptors to reads}
\description{Perform pairwise alignments of adaptor sequences to the ends of the read sequences.}

\usage{
adaptorAlign(adaptor1, adaptor2, reads, tolerance=100, 
    gapOpening=1, gapExtension=5, match=1, mismatch=0) 
}

\arguments{
\item{adaptor1, adaptor2}{A string or DNAString object containing the 5'-to-3' sequences of the adaptors on each end of the read.}
\item{reads}{A DNAStringSet object containing the read sequences.
A QualityScaledDNAStringSet can also be supplied containing quality scores.}
\item{tolerance}{An integer scalar specifying the region at the ends of the reads to search for the adaptors.}
\item{gapOpening, gapExtension}{Arguments to pass to \code{\link{pairwiseAlignment}}.}
\item{match, mismatch}{Arguments to pass to \code{\link{nucleotideSubstitutionMatrix}} for use in \code{\link{pairwiseAlignment}}.
These are ignored if \code{reads} contains quality scores.}
} 

\value{
A named list with the following components:
\describe{
\item{\code{reads}:}{A DNAStringSet or QualityScaledDNAStringSet object containing the read sequences.
If necessary, reverse-complementing is performed on specific entries such that the 5' end is always the one that aligns best to adaptor 1.}
\item{\code{adaptor1}:}{A DataFrame of alignment information for adaptor 1, including the pairwise alignment strings for the adaptor and read sequence; 
the start and end positions of the alignment on the read; the alignment score; and the qualities of the read subsequence, if they were originally supplied.
The read alignment string corresponds to the forward strand of the output \code{reads}, \emph{not} the input \code{reads}.
Alignment positions also refer to the forward strand of the output \code{reads}.
}
\item{\code{adaptor2}:}{A similar DataFrame of alignment information for adaptor 2.
The read alignment string corresponds to the \emph{reverse} strand of the output \code{reads}.
Alignment positions still refer to the forward strand of the output \code{reads}.
}
\item{\code{reversed}:}{A logical vector indicating whether each read has been reverse-complemented in the output.}
\item{\code{parameters}:}{A list containing the values of the arguments supplied to \code{adaptorAlign}, e.g., gap opening/extension penalties, \code{tolerance}.}
}
}

\details{
This function performs a local-global alignment of each adaptor to the ends of the read sequences.
This aims to identify the adaptors within the read sequence for trimming and/or UMI identification with \code{\link{chopReads}} and \code{\link{umiExtract}}.
Alignment is performed using the \code{\link{pairwiseAlignment}} function from the \pkg{Biostrings} package.
The default alignment parameters are chosen to account for the Nanopore's relatively high indel rate.

The adaptor sequences should be provided 5'-to-3', corresponding to the sequence on the ssDNA molecule.
Reverse-complementing is performed automatically within the function to search both strands of the read sequence.
There should usually be some experimental distinction between adaptor 1 and 2.
For example, in RNA sequencing, adaptor 1 may be used for reverse transcription while adaptor 2 may be used for second-strand synthesis.

Names are automatically generated for all reads if they are not supplied in the input \code{reads}.
These are stored in the names of the output \code{reads} as well as the row names of the output \code{adaptor1} and \code{adaptor2} DataFrames.
}

\section{Flipping reads by strand}{
This function ensures that the 5' end of the output read sequence is the one that aligns best to adaptor 1.
To this end, reads are reverse-complemented, quality strings are reversed and alignment positions are adjusted accordingly.
Thus, adaptor 1 will always be reported as aligning to the forward strand, while adaptor 2 will always be reported as aligning to the reverse strand.
This simplifies downstream analyses and comparions between read sequences, e.g., in \code{\link{consensusReadSeq}}.

The choice of the best orientation is performed by considering the sum of bounded alignment scores for each orientation.
To explain:
\enumerate{
\item Let us conceptually combine the two alignments of adaptor 1 and 2 to the (input) read sequence.
This can be treated as a single combined alignment, by imagining a stretch of \code{N}s between the two adaptors that contribute nothing to the score.
\item The score of the combined alignment is equal to the sum of scores for the individual alignments, due to the additive nature of the Smith-Waterman procedure.
However, if either individual score is negative, we set it to zero rather than adding a negative value.
This mimics the effect of local alignment, where an adaptor is ignored on either end if it does not improve the score.
\item We do the same for the scores for the reverse-complemented read sequence.
If this is higher than the original score, we decide that the read needs to be reverse-complemented.
}
The same approach is used to compute the scores that are returned by \code{\link{tuneAlignment}}.
}


\author{
Florian Bieberich,
Aaron Lun
}

\seealso{
\code{\link{chopReads}},
\code{\link{umiExtract}}
}

\examples{
a1 <- "AACGGGTCG"
a2 <- "CGTGCTGCATCG"
reads <- DNAStringSet(c(
    "AACGGGTCGCTCGACGACTACGTCCGATGCAGCACG",
    "CGTGCTGCATCGGACGTAGTCGTCGAGCGACCCGTT"))
(out <- adaptorAlign(adaptor1=a1, adaptor2=a2, reads=reads))
}
