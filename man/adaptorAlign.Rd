\name{adaptorAlign}
\alias{adaptorAlign}

\title{Align adaptors to reads}
\description{Perform pairwise alignments of adaptor sequences to the ends of the read sequences.}

\usage{
adaptorAlign(adaptor1, adaptor2, reads, quality = NULL, 
    gapOpening=1, gapExtension=5, match=2, mismatch=-5, 
    tolerance=100) 
}

\arguments{
\item{adaptor1, adaptor2}{A string or DNAString object containing the 5'-to-3' sequences of the adaptors on each end of the read.}
\item{reads}{A DNAStringSet object containing the read sequences.}
\item{quality}{A BStringSet or PhredQuality object containing Phred quality scores for each read.}
\item{gapOpening, gapExtension}{Arguments to pass to \code{\link{pairwiseAlignment}}.}
\item{match, mismatch}{Arguments to pass to \code{\link{nucleotideSubstitutionMatrix}} for use in \code{\link{pairwiseAlignment}}.}
\item{tolerance}{An integer scalar specifying the region at the ends of the reads to search for the adaptors.}
} 

\value{
A named list with the following components:
\describe{
\item{\code{reads}:}{A DNAStringSet object of the read sequences.
If necessary, reverse-complementing is performed on specific entries such that the 5' end is always the one that aligns best to adaptor 1.}
\item{\code{adaptor1}:}{A DataFrame of alignment information for adaptor 1, including the pairwise alignment strings for the adaptor (subject) and read sequence (pattern); 
and the start and end positions of the alignment on the adaptor and read.
The pattern alignment string corresponds to the forward strand of the output \code{reads}, \emph{not} the input \code{reads}.}
\item{\code{adaptor2}:}{A similar DataFrame of alignment information for adaptor 2.
The pattern alignment string corresponds to the reverse strand of the output \code{reads}.}
\item{\code{quality}:}{An object of the same type as the input \code{quality}, with entries reversed as described for the output \code{reads}.}
\item{\code{reversed}:}{A logical vector indicating whether each read has been reverse-complemented in the output.}
}
}

\details{
This function performs a local-global alignment of each adaptor to the ends of the read sequences.
This aims to identify the adaptors within the read sequence for trimming and/or UMI identification with \code{\link{chopReads}} and \code{\link{umiExtract}}.
Alignment is performed using the \code{\link{pairwiseAlignment}} function from the \pkg{Biostrings} package.
The default parameters are chosen to account for the Nanopore's relatively high indel rate.

The adaptor sequences should be provided 5'-to-3', corresponding to the sequence on the ssDNA molecule.
Reverse-complementing is performed automatically within the function to search both strands of the read sequence.
There may or may not be any distinction between adaptor 1 and 2, depending on whether the library preparation is strand-aware.
For example, in RNA sequencing, adaptor 1 may be used for reverse transcription while adaptor 2 may be used for second-strand synthesis.

This function ensures that the 5' end of the output read sequence is the one that aligns best to adaptor 1.
To this end, reads are reverse-complemented, quality strings are reversed and alignment positions are adjusted accordingly.
Thus, adaptor 1 will always be reported as aligning to the forward strand, while adaptor 2 will always be reported as aligning to the reverse strand.
This simplifies downstream analyses when the protocol is strand-aware (and has no effect when the protocol is not).

Names are automatically generated for all reads if they are not supplied in the input \code{reads}.
These are stored in the names of the output \code{reads} as well as the row names of the output \code{adaptor1} and \code{adaptor2} DataFrames.
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
