# This script tests various aspects of the consensus sequence machinery.
# library(sarlacc); library(testthat); source("test-consensus.R")

BASICFUN <- function(current, min.coverage, pseudo.count) {
    nucleotides <- c("A", "C", "G", "T")
    x <- consensusMatrix(current) 
    x <- x[c(nucleotides, "N"),,drop=FALSE] 
    base.exists <- colSums(x) >= min.coverage * length(current)
    x <- x[nucleotides,base.exists,drop=FALSE]
    
    # Defining the consensus sequence based on relative majority coverage
    # (at least one base must be above 25%, as we only keep 'x' in "ATCG")
    x <- t(x)
    majority.base <- max.col(x, ties.method="first")
    conseq <- nucleotides[majority.base]
    consensus <- paste(conseq, collapse="")
    
    # Computing the Phred score.
    n.chosen <- numeric(length(conseq))
    for (j in nucleotides) { 
        chosen <- conseq==j
        n.chosen[chosen] <- x[chosen,j]
    }
    phred <- 1 - (n.chosen + pseudo.count/4)/(rowSums(x)+pseudo.count)
    return(list(consensus, phred))
}

test_that("basic consensus formation works properly", {
    alignments <- c("AAAAGAAAAA-AAATAAAA",
                    "ACACA-AAAA--AAT-AGA",
                    "GA-AG-C-A-T-AAT-AAA",
                    "AT-AG-T-AGTAAGA-AGA",
                    "-AAAGAT-AGTCAGA-AGA",
                    "AGAAAAT-AGAAATA-AGA") # variety of errors; indels and/or substitutions

    library(Biostrings)
    current <- DNAStringSet(alignments)
    expect_identical(BASICFUN(current, min.coverage=0.6, pseudo.count=1),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.6, 1))
    expect_identical(BASICFUN(current, min.coverage=0.6, pseudo.count=2),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.6, 2))
    expect_identical(BASICFUN(current, min.coverage=0.2, pseudo.count=2),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.2, 2))
    expect_identical(BASICFUN(current, min.coverage=0.9, pseudo.count=1),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.9, 1))

    # Throwing in variable numbers of Ns, with and without deletions.
    # (note that N's default to A's, but obviously at very low quality).
    alignments <- c("NAAAAANNN",
                    "NNAANA---",
                    "NNNANNN--",
                    "NNNN--NN-")
    current <- DNAStringSet(alignments)
    expect_identical(BASICFUN(current, min.coverage=0.6, pseudo.count=1),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.6, 1))
    expect_identical(BASICFUN(current, min.coverage=0.6, pseudo.count=2),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.6, 2))
    expect_identical(BASICFUN(current, min.coverage=0.2, pseudo.count=2),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.2, 2))
    expect_identical(BASICFUN(current, min.coverage=0.9, pseudo.count=1),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.9, 1))
})
