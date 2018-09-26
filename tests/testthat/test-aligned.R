# This tests the C++ function to extract the aligned subject and pattern.
# library(sarlacc); library(testthat); source("test-aligned.R")

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
adaptor <- DNAString("AAAAGGGGCCCCTTTT")

##########################################

test_that("alignment function is behaving correctly", {
    aln <- pairwiseAlignment(pattern=reads, subject=adaptor, type="local-global", gapOpening=1)
    obs <- sarlacc:::.bplalign(reads=reads, adaptor=adaptor, type="local-global", gapOpening=1)

    expect_identical(obs$score, score(aln))
    expect_identical(obs$read, as.character(alignedPattern(aln)))
    expect_identical(obs$adaptor, as.character(alignedSubject(aln)))

    true.pos <- gregexpr("[AGCT]", read.seq)
    true.start <- unlist(lapply(true.pos, min))
    true.end <- unlist(lapply(true.pos, max))
    expect_identical(obs$start, true.start)
    expect_identical(obs$end, true.end)
})

test_that("alignment function behaves with quality strings", {
    quals <- gsub("[acgt]","!", gsub("[ACGT]", "7", read.seq))
    qreads <- QualityScaledDNAStringSet(reads, PhredQuality(quals))
    obs <- sarlacc:::.bplalign(reads=qreads, adaptor=adaptor, type="local-global", gapOpening=1)

    expect_identical(nchar(obs$quality), nchar(gsub("-", "", obs$read)))
    expect_identical(nchar(obs$quality), obs$end - obs$start + 1L)
    expect_true(all(grepl("^7.*7$", obs$quality)))
})

test_that("alignment behaves correctly around ambiguous bases", {
    adaptorN <- DNAString("AAAAGGNNNNCCTTTT")

    # No quality strings.
    args <- sarlacc:::.setup_alignment_args(has.quality=FALSE, gapOpening=5, gapExtension=1, match=5, mismatch=0)
    expect_true(!is.null(args$substitutionMatrix))
    expect_true(is.null(args$fuzzyMatrix))

    ref <- do.call(sarlacc:::.bplalign, c(list(reads=reads, adaptor=adaptor), args))
    out <- do.call(sarlacc:::.bplalign, c(list(reads=reads, adaptor=adaptorN), args))
    expect_identical(ref$read, out$read)
    expect_identical(ref$start, out$start)
    expect_identical(ref$end, out$end)

    # With qualities.
    quals <- PhredQuality(strrep("9", width(reads)))
    qreads <- QualityScaledDNAStringSet(reads, quals)
    qargs <- sarlacc:::.setup_alignment_args(has.quality=TRUE, gapOpening=5, gapExtension=1)
    expect_true(is.null(qargs$substitutionMatrix))
    expect_true(!is.null(qargs$fuzzyMatrix))

    qref <- do.call(sarlacc:::.bplalign, c(list(reads=qreads, adaptor=adaptor), qargs))
    qout <- do.call(sarlacc:::.bplalign, c(list(reads=qreads, adaptor=adaptorN), qargs))
    expect_identical(ref$read, qout$read)
    expect_identical(ref$start, qout$start)
    expect_identical(ref$end, qout$end)

    # For comparison, to show that the 'fuzzyMatrix' is necessary:
    bad.args <- qargs
    bad.args$fuzzyMatrix <- NULL
    bad.out <- do.call(sarlacc:::.bplalign, c(list(reads=qreads, adaptor=adaptorN), bad.args))
    expect_false(identical(bad.out$read, qref$read))
    expect_true(all(bad.out$score <= qref$score + 1e-8))
})

test_that("alignment function handles multiple cores and empty inputs", {
    # Multiple cores:
    ref <- sarlacc:::.bplalign(adaptor=adaptor, reads=reads)
    multi2 <- sarlacc:::.bplalign(adaptor=adaptor, reads=reads, BPPARAM=MulticoreParam(2))
    expect_identical(multi2, ref)
    multi3 <- sarlacc:::.bplalign(adaptor=adaptor, reads=reads, BPPARAM=MulticoreParam(3))
    expect_identical(multi3, ref)

    # Empty inputs:
    empty <- sarlacc:::.bplalign(adaptor=adaptor, reads=reads[0])
    expect_identical(ref[0,], empty)
    empty <- sarlacc:::.bplalign(adaptor=adaptor, reads=reads[0], scoreOnly=TRUE)
    expect_identical(empty, numeric(0))
})

##########################################

test_that("quality assigner works correctly", {
    # As character strings:
    Q <- sarlacc:::.assign_qualities(read.seq)
    expect_s4_class(Q, "QualityScaledDNAStringSet")
    expect_identical(seq_along(Q), grep("^~+$", as.character(quality(Q))))
    expect_identical(toupper(read.seq), as.character(Q))

    # As a DNAStringSet.
    Q2 <- sarlacc:::.assign_qualities(reads)
    expect_s4_class(Q2, "QualityScaledDNAStringSet")
    expect_identical(quality(Q), quality(Q2))

    # As a QSDS:
    quals <- gsub("[acgt]","!", gsub("[ACGT]", "7", read.seq))
    qreads <- QualityScaledDNAStringSet(reads, PhredQuality(quals))
    Q3 <- sarlacc:::.assign_qualities(qreads)
    expect_identical(Q3, qreads)
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

##########################################

test_that("overall adaptorAlign function works correctly", {
    myread <- "AAAAAAACGTACGTACGTGGGGGGG"
    revread <- as.character(reverseComplement(DNAString(myread)))

    # Checking that the flipping of reads works correctly:
    out <- adaptorAlign("AAAAAAA", "CCCCCCC", c(myread, revread))
    expect_identical(out$reversed, c(FALSE, TRUE))
    expect_identical(out$reads[1], DNAStringSet(myread))
    expect_identical(out$reads[1], out$reads[2])

    expect_identical(out$adaptor1$start[1], 1L)
    expect_identical(out$adaptor1$end[1], 7L)
    expect_identical(out$adaptor2$start[1], nchar(myread))
    expect_identical(out$adaptor2$end[1], nchar(myread)-7L+1L)

    expect_identical(out$adaptor1[1,], out$adaptor1[2,])
    expect_identical(out$adaptor2[1,], out$adaptor2[2,])

    # Checking that quality reporting works correctly:
    myqual <- "1234567890ABCDEFGHIJKLMNO"
    QSDS <- QualityScaledDNAStringSet(c(myread, revread), PhredQuality(rep(myqual, 2)))
    out <- adaptorAlign("AAAAAAA", "CCCCCCC", QSDS)
    expect_identical(out$reversed, c(FALSE, TRUE))
    expect_identical(out$reads[1], QSDS[1])
    expect_identical(out$reads[2], reverseComplement(QSDS[2]))

    expect_identical(as.character(out$adaptor1$quality[1]), substr(myqual, 1, 7))
    expect_identical(as.character(out$adaptor1$quality[2]), reverse(substr(myqual, start=nchar(myqual) - 7 + 1, stop=nchar(myqual))))
    expect_identical(as.character(out$adaptor2$quality[1]), as.character(out$adaptor1$quality[2]))
    expect_identical(as.character(out$adaptor2$quality[2]), as.character(out$adaptor1$quality[1]))

    # Handles empty inputs.
    empty <- adaptorAlign("AAAAAAA", "CCCCCCC", QSDS[0])
    expect_identical(nrow(empty), 0L)
})
