# This tests the C++ functions for pairwise sequence alignment.
# library(sarlacc); library(testthat); source("test-aligned.R")

library(Biostrings)
adaptor <- "AAAAGGGGCCCCTTTT"

read.seq <- c("AAAAGGGGCCCCTTTT", # identical
              "acgtacgtacgtAAAAGGGGCCCCTTTT", # insertion at the start
              "AAAAGGGGCCCCTTTTacgtacgtacgt", # insertion at the end
              "GGGGCCCCTTTT", # deletion at the start
              "AAAAGGGGCCCC", # deletion at the end

              "acgtacgtacgtAAAAGGGGCCCCTTTTacgtacgtacgt", # insertion at the start, insertion at the end
              "acgtacgtacgtAAAAGGGGCCCC", # insertion at the start, deletion at the end
              "GGGGCCCCTTTTacgtacgtacgt", # deletion at the start, insertion at the end
              "GGGGCCCC", # deletion at the start, deletion at the end

              "AAAAGGGGacgtCCCCTTTT", # insertion in the middle
              "AAAAGGCCTTTT") # deletion in the middle

reads <- DNAStringSet(read.seq)

set.seed(1000)
quals <- do.call(c, lapply(read.seq, FUN=function(x) { PhredQuality(rbeta(nchar(x), 1, 10)) }))
qreads <- QualityScaledDNAStringSet(reads, quals)
   
REF_ALIGN <- function(R, A, gapOpening, gapExtension, type="local-global") {
    pairwiseAlignment(pattern=R, 
        subject=QualityScaledDNAStringSet(DNAString(A), PhredQuality(rep(0, nchar(A)))),
        fuzzyMatrix=nucleotideSubstitutionMatrix(), 
        type=type,
        gapExtension=gapExtension, gapOpening=gapOpening)            
}

test_that("alignment scores and positions are computed correctly", {
    ref <- REF_ALIGN(qreads, adaptor, gapOpening=5, gapExtension=1)

    # Alignment scores match up (mostly, as 'error' for the reference is non-zero).
    out <- sarlacc:::.align_and_extract(adaptor, qreads, gap.opening=5, gap.extension=1, subseq.starts=integer(0), subseq.ends=integer(0))
    expect_equal(score(ref), out$score, tol=0.0001) 

    # Alignment positions match up.
    read0 <- pattern(ref)
    expect_identical(start(read0), out$start)
    expect_identical(end(read0), out$end)

    # Function behaves with empty adaptor.
    out <- sarlacc:::.align_and_extract("", qreads, gap.opening=5, gap.extension=1, subseq.starts=integer(0), subseq.ends=integer(0))
    expect_identical(out$score, numeric(nrow(out)))
    expect_identical(out$start, integer(nrow(out)))
    expect_identical(out$end, integer(nrow(out)))

    out <- sarlacc:::.align_and_extract(adaptor, subseq(qreads, start=1, width=0), gap.opening=5, gap.extension=1, subseq.starts=integer(0), subseq.ends=integer(0))
    expect_identical(out$score, rep(-nchar(adaptor) - 5, nrow(out)))
    expect_identical(out$start, integer(nrow(out)))
    expect_identical(out$end, integer(nrow(out)))
})

test_that("affine gap penalties are handled correctly", {
    # Affine gap penalties complicate the DP matrix, which needs to consider
    # whether a gap opening in the previous position (which would result in a
    # suboptimal score there) might lead to an optimal score in the current position.

    # Here, we set up a scenario where one mismatch penalty is less damaging than
    # a single gap but multiple mismatches are more damaging than a gap of the same length.
    # We first do so with gaps in the read versus the reference.
    a1 <- "AAACCCAAATTTAAA"
    qread <- QualityScaledDNAStringSet("AAAAAAAAA", PhredQuality(strrep("+", 9))) 
    ref <- REF_ALIGN(qread, a1, 5, 1)

    out <- sarlacc:::.align_and_extract(a1, qread, gap.opening=5, gap.extension=1, subseq.starts=integer(0), subseq.ends=integer(0))
    expect_equal(score(ref), out$score, tol=0.0001)
    expect_identical(start(pattern(ref)), out$start)
    expect_identical(end(pattern(ref)), out$end)

    # We repeat this process with gaps in the reference versus the read.
    a1 <- "AAAAAA"
    qread <- QualityScaledDNAStringSet("AAACCCAAA", PhredQuality(strrep("+", 9))) 
    ref <- REF_ALIGN(qread, a1, 5, 1)

    out <- sarlacc:::.align_and_extract(a1, qread, gap.opening=5, gap.extension=1, subseq.starts=integer(0), subseq.ends=integer(0))
    expect_equal(score(ref), out$score, tol=0.0001)
    expect_identical(start(pattern(ref)), out$start)
    expect_identical(end(pattern(ref)), out$end)
})

test_that("alignment extraction works correctly", {
    ref <- REF_ALIGN(qreads, adaptor, gapOpening=5, gapExtension=1)
    refR <- as.character(alignedPattern(ref))
    refA <- as.character(alignedSubject(ref))
    gapsA <- lapply(strsplit(refA, ""), FUN=function(x) cumsum(x!="-"))

    # Extracts the correct components.
    possibilities <- combn(nchar(adaptor)+1L, 2)
    possibilities[2,] <- possibilities[2,] - 1L 
    possibilities <- possibilities[,sample(ncol(possibilities))]

    out <- sarlacc:::.align_and_extract(adaptor, qreads, gap.opening=5, gap.extension=1, subseq.starts=possibilities[1,], subseq.ends=possibilities[2,])

    # Comparing to a reference implementation.
    for (i in seq_len(ncol(possibilities))) {
        observed <- as.character(out$subseq[,i])
        curstart <- possibilities[1,i]
        curend <- possibilities[2,i]

        collected <- character(length(observed))
        for (j in seq_along(gapsA)) {
            collected[j] <- substr(refR[j], min(which(gapsA[[j]]==curstart)), min(which(gapsA[[j]]==curend)))
        }
        collected <- gsub("-", "", collected)
        expect_identical(observed, collected)
    }
})

test_that("subsequence finder behaves correctly around ambiguous bases", {
    expect_identical(sarlacc:::.setup_subseqs("AAAAGGNNNNCCTTTT"), list(starts=7L, ends=10L))
    expect_identical(sarlacc:::.setup_subseqs("AAAAGGYYYYCCTTTT"), list(starts=7L, ends=10L))
    expect_identical(sarlacc:::.setup_subseqs("AAAAGGCCTTTTRRRR"), list(starts=13L, ends=16L))
})

test_that("sequence front/back getter works correctly", {
    extremes <- sarlacc:::.get_front_and_back(reads, tolerance=10)
    expect_identical(as.character(extremes$front), toupper(substr(read.seq, 1, 10)))
    expect_identical(as.character(reverseComplement(extremes$back)), toupper(substr(read.seq, nchar(read.seq) - 10 + 1, nchar(read.seq))))

    # Handles past-the range.
    extremes <- sarlacc:::.get_front_and_back(reads, tolerance=10000)
    expect_identical(as.character(extremes$front), toupper(read.seq))
    expect_identical(as.character(reverseComplement(extremes$back)), toupper(read.seq))
})

test_that("overall adaptorAlign function works correctly", {
    myread <- "AACGTAACGTACGTACGTGGGGGGG"
    myqual <- "1234567890ABCDEFGHIJKLMNO"
    QSDS <- QualityScaledDNAStringSet(myread, PhredQuality(myqual))
    QSDS <- c(QSDS, reverseComplement(QSDS))

    tmp <- tempfile(fileext=".fastq")
    names(QSDS) <- c("X", "Y")
    writeXStringSet(QSDS, qualities=quality(QSDS), format="fastq", filepath=tmp)

    # Checking that the flipping of reads works correctly:
    out <- adaptorAlign("AANNNAA", "CCCCCCC", tmp)

    expect_identical(out$reversed, c(FALSE, TRUE))
    expect_identical(nrow(unique(out$adaptor1)), 1L) # otherwise identical.
    expect_identical(nrow(unique(out$adaptor2)), 1L) # otherwise identical.

    expect_identical(out$adaptor1$start[1], 1L)
    expect_identical(out$adaptor1$end[1], 7L)
    expect_identical(out$adaptor2$start[1], nchar(myread))
    expect_identical(out$adaptor2$end[1], nchar(myread)-7L+1L)

    # Handles empty inputs.
    writeXStringSet(QSDS[0], qualities=quality(QSDS)[0], format="fastq", filepath=tmp)
    empty <- adaptorAlign("AAAAAAA", "CCCCCCC", tmp)
    expect_identical(nrow(empty), 0L)
})

test_that("global alignment scores are computed correctly", {
    ref <- REF_ALIGN(qreads, adaptor, gapOpening=5, gapExtension=1, type="global")

    # Alignment scores match up (mostly, as 'error' for the reference is non-zero).
    out <- sarlacc:::.align_BA_internal(qreads, adaptor, gap.opening=5, gap.extension=1)
    expect_equal(score(ref), out, tol=0.0001) 

    # Function behaves with empty adaptor.
    out <- sarlacc:::.align_BA_internal(qreads, "", gap.opening=5, gap.extension=1)
    expect_identical(out, - width(qreads)  - 5)

    out <- sarlacc:::.align_BA_internal(subseq(qreads, start=1, width=0), adaptor, gap.opening=5, gap.extension=1)
    expect_identical(out, rep(-nchar(adaptor) - 5, length(out)))
})
