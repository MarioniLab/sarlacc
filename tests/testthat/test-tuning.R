# This tests the tuneAlignment() function.
# library(testthat); library(sarlacc); source("test-tuning.R")

a1 <- "CGTACGACGAT"
a2 <- "TCGAGCGTTAC"
reads <- c(
    "CGTACGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGCTCGA",
    "CGTCGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGCTCGA",
    "CGTACGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGCCGA",
    "CGTACGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGCTCGA",
    "CGCTACGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGGTCGA"
)

test_that("tuneAlignments works correctly without quality scores", {
    out <- tuneAlignment(a1, a2, reads, gapOp.range=c(4, 5), gapExt.range=c(1,2))
    expect_true(!is.na(out$parameter$match))
    expect_true(!is.na(out$parameter$mismatch))

    args <- list(substitutionMatrix=nucleotideSubstitutionMatrix(out$parameter$match, out$parameter$mismatch),
        gapOpening=out$parameter$gapOpening, gapExtension=out$parameter$gapExtension, type="local-global")

    aln1F <- do.call(pairwiseAlignment, c(list(pattern=DNAStringSet(reads), subject=DNAString(a1)), args))
    aln2F <- do.call(pairwiseAlignment, c(list(pattern=DNAStringSet(reads), subject=reverseComplement(DNAString(a2))), args))
    aln1R <- do.call(pairwiseAlignment, c(list(pattern=DNAStringSet(reads), subject=reverseComplement(DNAString(a1))), args))
    aln2R <- do.call(pairwiseAlignment, c(list(pattern=DNAStringSet(reads), subject=DNAString(a2)), args))

    Fscores <- pmax(score(aln1F), 0) + pmax(score(aln2F), 0)
    Rscores <- pmax(score(aln1R), 0) + pmax(score(aln2R), 0)
    expect_equal(out$scores$reads, pmax(Fscores, Rscores))

    expect_true(min(out$scores$reads) > max(out$scores$scrambled))
})


test_that("tuneAlignments works correctly with quality scores", {
    qreads <- QualityScaledDNAStringSet(DNAStringSet(reads),
        PhredQuality(strrep("2", nchar(reads))))
    qa1 <- QualityScaledDNAStringSet(DNAStringSet(a1), PhredQuality(strrep("~", nchar(a1))))
    qa2 <- QualityScaledDNAStringSet(DNAStringSet(a2), PhredQuality(strrep("~", nchar(a2))))

    out <- tuneAlignment(qa1, qa2, qreads, gapOp.range=c(4, 5), gapExt.range=c(1,2))
    expect_true(is.na(out$parameter$match))
    expect_true(is.na(out$parameter$mismatch))
                
    args <- list(fuzzyMatrix=nucleotideSubstitutionMatrix(), gapOpening=out$parameter$gapOpening, gapExtension=out$parameter$gapExtension, type="local-global")

    aln1F <- do.call(pairwiseAlignment, c(list(pattern=qreads, subject=qa1), args))
    aln2F <- do.call(pairwiseAlignment, c(list(pattern=qreads, subject=reverseComplement(qa2)), args))
    aln1R <- do.call(pairwiseAlignment, c(list(pattern=qreads, subject=reverseComplement(qa1)), args))
    aln2R <- do.call(pairwiseAlignment, c(list(pattern=qreads, subject=qa2), args))

    Fscores <- pmax(score(aln1F), 0) + pmax(score(aln2F), 0)
    Rscores <- pmax(score(aln1R), 0) + pmax(score(aln2R), 0)
    expect_equal(out$scores$reads, pmax(Fscores, Rscores))

    expect_true(min(out$scores$reads) > max(out$scores$scrambled))
})

test_that("tuneAlignments fails gracefully with no reads", {
    blah <- tuneAlignment(a1, a2, reads[0])
    expect_identical(rep(NA_integer_, 4), unname(unlist(blah$parameters)))
    expect_identical(integer(2), unname(lengths(blah$scores)))
})

test_that("tied FDR calculator works as expected", {
    expect_equal(1, sarlacc:::.tied_overlap(1:10, 1:10 - 10))
    expect_equal(0.5, sarlacc:::.tied_overlap(1:10, 1:10))
    expect_equal(0.55, sarlacc:::.tied_overlap(1:10, 1:10 - 0.5))
    expect_equal(0.45, sarlacc:::.tied_overlap(1:10, 1:10 + 0.5))
    expect_equal(0, sarlacc:::.tied_overlap(1:10, 1:10 + 10))
})
