# This tests the C++ functions for pairwise sequence alignment.
# library(sarlacc); library(testthat); source("test-general-align.R")

library(Biostrings)
   
nucleotides <- c("A", "C", "G", "T")
GENERATOR <- function() {
    ref <- paste(sample(nucleotides, 50, replace=TRUE), collapse="")
    obs <- paste(sample(nucleotides, sample(20:80, 1), replace=TRUE), collapse="")
    
    quals <- PhredQuality(10^-runif(nchar(obs), 1, 5))
    q <- QualityScaledDNAStringSet(obs, quals)
    r <- QualityScaledDNAStringSet(ref, PhredQuality(strrep("~", nchar(ref))))
    list(ref=r, query=q)
}

STRIPPER <- function(qalign1, ralign1, qalign2, ralign2, quals) 
# We need a special comparison function; alignment strings
# are not unique as there are multiple best backtrack paths.
{
    expect_identical(nchar(ralign1), nchar(ralign2))
    expect_identical(nchar(qalign1), nchar(qalign2))

    # Check that the underlying sequences are the same.
    R1 <- gsub("-", "", ralign1)
    Q1 <- gsub("-", "", qalign1)
    R2 <- gsub("-", "", ralign2)
    Q2 <- gsub("-", "", qalign2)
    expect_identical(R1, R2)
    expect_identical(Q1, Q2)

    # Check that the edit distances are the same.
    rvec1 <- strsplit(ralign1, "")[[1]]
    qvec1 <- strsplit(qalign1, "")[[1]]
    rvec2 <- strsplit(ralign2, "")[[1]]
    qvec2 <- strsplit(qalign2, "")[[1]]
    expect_identical(sum(rvec1!=qvec1), sum(rvec2!=qvec2))

    # Check that it's the same query bases that are mismatched;
    # or if they differ, then it's those bases have the same quality
    # (which would also result in equally best backtrack paths).
    mm1 <- qvec1!=rvec1 & rvec1!='-' & qvec1!='-'
    mm2 <- qvec2!=rvec2 & rvec2!='-' & qvec2!='-'
    mismatched1 <- cumsum(qvec1!='-')[mm1]
    mismatched2 <- cumsum(qvec2!='-')[mm2]
    expect_identical(length(mismatched1), length(mismatched2))

    notsame <- mismatched1!=mismatched2
    quals <- strsplit(quals, "")[[1]]
    expect_identical(quals[mismatched1[notsame]], quals[mismatched2[notsame]])

    NULL
}

set.seed(32000)
test_that("generalAlign works correctly for global alignments", {
    for (i in seq_len(50)) {
        sim <- GENERATOR()
        r <- sim$ref
        q <- sim$query
        quals <- as.character(quality(q))

        aln <- pairwiseAlignment(q, r, gapOpening=5, gapExtension=1)
        out <- qualityAlign(q, r, gapOpening=5)
        expect_equal(score(aln), out$score, tol=1e-5)
        STRIPPER(out$query, out$reference, as.character(alignedPattern(aln)), as.character(alignedSubject(aln)), quals)

        aln <- pairwiseAlignment(q, r, gapOpening=20, gapExtension=1)
        out <- qualityAlign(q, r, gapOpening=20)
        expect_equal(score(aln), out$score, tol=1e-5)
        STRIPPER(out$query, out$reference, as.character(alignedPattern(aln)), as.character(alignedSubject(aln)), quals)

        aln <- pairwiseAlignment(q, r, gapOpening=1, gapExtension=1)
        out <- qualityAlign(q, r, gapOpening=1)
        expect_equal(score(aln), out$score, tol=1e-5)
        STRIPPER(out$query, out$reference, as.character(alignedPattern(aln)), as.character(alignedSubject(aln)), quals)
    }
})

set.seed(32001)
test_that("generalAlign works correctly for multiple queries", {
    rsim <- GENERATOR()
    ref <- rsim$ref

    collected <- results <- vector("list", 20)
    for (i in seq_len(20)) {
        collected[[i]] <- GENERATOR()$query
        results[[i]] <- qualityAlign(collected[[i]], ref)
    }

    output <- qualityAlign(do.call(c, collected), ref)
    expect_identical(as.data.frame(output), as.data.frame(do.call(rbind, results)))
})

