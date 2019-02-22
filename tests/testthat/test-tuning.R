# This tests the tuneAlignment() function.
# library(testthat); library(sarlacc); source("test-tuning.R")

a1 <- "CGTACGACGAT"
a2 <- "TCGAGCGTTAC"
reads <- c(
    "CGTACGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGCTCGA",
    "CGTCGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGCTCGA",
    "CGTACGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGCCGA",
    "CGTACGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGCTCGA",
    "CGCTACGACGATGACTGATCGATCGTAGTTCATCGACGATGTAACGGTCGA",
    "GTACGACGATTTACTGATCGATCGTAGTTCATCGACGATGTAACGCTCGA",
    "CGTCGACGATGACTGGGCGATCGTAGTTCATCGACGATGTAACGCTCGA",
    "TACGACGATGACTGATCATCGTAAAAATTCATCGATGTAACGCCGA",
    "CGTACGACGATGACTGATCGATCGTAGTTCACCCCGATGTAACGCTCGA",
    "CGCTACGACGATGACTGCGATCGTAGTTCATAAAAATGTAACGGTCGA"
)

qreads <- QualityScaledDNAStringSet(DNAStringSet(reads),
    PhredQuality(strrep("~", nchar(reads))))

tmp <- tempfile(fileext=".fastq")
names(qreads) <- sprintf("READ_%i", seq_along(reads))
writeXStringSet(qreads, qualities=quality(qreads), format="fastq", filepath=tmp)

test_that("tuneAlignments works correctly with quality scores", {
    out <- tuneAlignment(a1, a2, tmp, gapOp.range=c(4, 5), gapExt.range=c(1,2))
    args <- list(fuzzyMatrix=nucleotideSubstitutionMatrix(), gapOpening=out$parameter$gapOpening, gapExtension=out$parameter$gapExtension, type="local-global")

    qa1 <- QualityScaledDNAStringSet(DNAStringSet(a1), PhredQuality(strrep("~", nchar(a1))))
    qa2 <- QualityScaledDNAStringSet(DNAStringSet(a2), PhredQuality(strrep("~", nchar(a2))))
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
    tmp2 <- tempfile(fileext=".fastq")
    writeXStringSet(qreads[0], qualities=quality(qreads[0]), format="fastq", filepath=tmp2)

    blah <- tuneAlignment(a1, a2, tmp2)
    expect_identical(rep(NA_integer_, 2), unname(unlist(blah$parameters)))
    expect_identical(integer(2), unname(lengths(blah$scores)))
})

test_that("tied FDR calculator works as expected", {
    expect_equal(1, sarlacc:::.tied_overlap(1:10, 1:10 - 10))
    expect_equal(0.5, sarlacc:::.tied_overlap(1:10, 1:10))
    expect_equal(0.55, sarlacc:::.tied_overlap(1:10, 1:10 - 0.5))
    expect_equal(0.45, sarlacc:::.tied_overlap(1:10, 1:10 + 0.5))
    expect_equal(0, sarlacc:::.tied_overlap(1:10, 1:10 + 10))
})
