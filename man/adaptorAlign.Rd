\name{adaptorAlign}
\alias{adaptorAlign}

\title{Align adaptors to reads}
\description{Perform pairwise alignments of adaptor sequences to the ends of the read sequences.}

\usage{
adaptorAlign(adaptor1, adaptor2, filepath, tolerance=100, gapOpening=1, 
    gapExtension=5, qual.type=c("phred", "solexa", "illumina"), 
    number=1e5, BPPARAM=SerialParam())
}

\arguments{
\item{adaptor1, adaptor2}{A string or DNAString object containing the 5'-to-3' sequences of the adaptors on each end of the read.}
\item{filepath}{A string containing the path to the FASTQ file, or a connection object to a FASTQ file.}
\item{tolerance}{An integer scalar specifying the length of the ends of the reads to search for adaptors.}
\item{gapOpening}{A numeric scalar specifying the gap opening penalty.}
\item{gapExtension}{A numeric scalar specifying the gap extension penalty.}
\item{qual.type}{String specifying the type of quality scores in \code{filepath}.}
\item{number}{Integer scalar specifying the number of records to read at once from the FASTQ file, see \code{?\link{FastqStreamer}}.}
\item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying whether alignment should be parallelized.}
} 

\value{
A \linkS4class{DataFrame} where each row corresponds to a read in \code{filepath} and is named with the read name.
It contains the following fields:
\describe{
\item{read.width}{Integer vector containing the length of each read.}
\item{\code{adaptor1}:}{A DataFrame of alignment information for adaptor 1, including the start and end positions of the alignment on the read; the alignment score; 
and the read subsequences aligned to any ambiguous bases in the adaptor.}
\item{\code{adaptor2}:}{A similar DataFrame of alignment information for adaptor 2.}
\item{\code{reversed}:}{A logical vector indicating whether each read has been reverse-complemented, see below.}
}

The \code{filepath}, \code{qual.type} and \code{tolerance} used are stored in the metadata of the overall DataFrame.
Alignment parameters and reference adaptor sequences are stored in the metadata of the each internal \code{adaptor*} DataFrame.
}

\details{
This function performs a local-global alignment of each adaptor to the ends of the read sequences (local for the read, global for the adaptor).
This aims to identify the adaptors within the read sequence for trimming and/or UMI identification.
Alignment is performed using a quality-aware algorithm akin to that in the \code{\link{pairwiseAlignment}} function from the \pkg{Biostrings} package.
The default alignment parameters are chosen to account for the Nanopore's relatively high indel rate.

The adaptor sequences should be provided 5'-to-3', corresponding to the sequence on the ssDNA molecule.
Reverse-complementing is performed automatically within the function to search both strands of the read sequence.
There should usually be some experimental distinction between adaptor 1 and 2.
For example, in RNA sequencing, adaptor 1 may be used for reverse transcription while adaptor 2 may be used for second-strand synthesis.

Any stretch of ambiguous IUPAC bases in the adaptor sequence is identified as a subsequence of interest.
This can represent unique molecular identifiers or barcode sequences that may be read-specific.
For each read, \code{adaptorAlign} will report the corresponding read subsequence that aligns to this stretch.
These are reported in the same orientation as the adaptor sequence and can be directly compared across reads.
}

\section{Flipping reads by strand}{
We standardize the alignment statistics in terms of the \dQuote{canonical cDNA orientation}.
In this canonical orientation, the sequence of adaptor 1 is located on the 5' end while the reverse complement of adaptor 2 is located on the 3' end.
Any reads observed in the opposite orientation are conceptually reverse-complemented, as marked by the \code{reversed} field in the output.
This standardization simplifies downstream analyses and comparisons between read sequences, e.g., in \code{\link{consensusReadSeq}}.

We decide whether or not a read is in the canonical orientation by considering the sum of bounded alignment scores for each orientation.
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

Start and end positions of the adaptor alignments are reported relative to the canonical orientation. 
\itemize{
\item For adaptor 1, all coordinates should generally be small, and the start position should be smaller than the end position.
This is because the start of adaptor 1 should align to the start of the read, while the end of the adaptor should align closer to the centre of the read.
\item For adaptor 2, all coordinates should generally be close to the read length, and the start position should be \emph{larger} than the end position.
This is because the start of adaptor 2 should align to the end of the read, while the end of the adaptor should align closer to the centre of the read.
}

The interpretation of the alignment strings is slightly different, as they are reported relative to the supplied adaptor sequences. 
The starts of the alignment strings always refer to the 5' end of the reference adaptor sequence, regardless of which adaptor is involved.
This simplifies downstream processing steps that extract the read sequence corresponding to a stretch of the adaptor sequence, e.g., \code{\link{umiExtract}}.
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
# Mocking up a data set.
a1 <- "AACGGGTCGNNNNNNNACGTACGTNNNNACGA" 
a2 <- "CGTGCTGCATCG"
fout <- tempfile(fileext=".fastq")
ref <- sarlacc:::mockReads(a1, a2, fout, nmolecules=3)

# Aligning it.
(out <- adaptorAlign(adaptor1=a1, adaptor2=a2, filepath=fout))
}
