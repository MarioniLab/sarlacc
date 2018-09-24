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
    adaptorN <- "AAAAGGNNNNCCTTTT"

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

    qref <- do.call(sarlacc:::.bplalign, c(list(reads=qreads, adaptor=adaptor), args))
    qout <- do.call(sarlacc:::.bplalign, c(list(reads=qreads, adaptor=adaptorN), args))
    expect_identical(ref$read, qout$read)
    expect_identical(ref$start, qout$start)
    expect_identical(ref$end, qout$end)
})

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
