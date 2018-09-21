\name{qualityMask}
\alias{qualityMask}

\title{Mask bases by quality}
\description{Mask out low-quality bases in a DNA sequence prior to downstream analyses.}

\usage{
qualityMask(seq, max.err) 
}

\arguments{
    \item{seq}{A \linkS4class{QualityScaledDNAStringSet} containing one or more sequences with quality strings.}
    \item{max.err}{Numeric scalar specifying the error probability threshold, above which bases are to be masked.}
}

\value{
A \linkS4class{DNAStringSet} containing masked DNA sequences.
All bases with error probabilities above \code{max.err} will be replaced with \code{N}s.
}

\author{
Aaron Lun
}

\examples{
X <- QualityScaledDNAStringSet(
    DNAStringSet(paste0("ACAGCTGCTGCAGCTGTG")),
    PhredQuality(rbeta(18, 1, 10))
)

qualityMask(X, max.err=0.5)
qualityMask(X, max.err=0.1)
qualityMask(X, max.err=0.01)
}

\seealso{
\code{\link{encoding}}
}
