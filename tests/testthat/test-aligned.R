# This tests the C++ function to extract the aligned subject and pattern.
# library(sarlacc); library(testthat); source("test-aligned.R")

test_that("alignment information extractor is behaving correctly", {
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
    
    aln <- pairwiseAlignment(pattern=DNAStringSet(read.seq), subject=DNAString("AAAAGGGGCCCCTTTT"), type="local-global", gapOpening=1)

    # Testing without quality strings.        
    out <- sarlacc:::.align_info_extractor(aln)
    expect_identical(out$read, as.character(alignedPattern(aln)))
    expect_identical(out$adaptor, as.character(alignedSubject(aln)))

    # Also checking what happens with quality strings.
    quals <- gsub("[acgt]","!", gsub("[ACGT]", "7", read.seq))
    out <- sarlacc:::.align_info_extractor(aln, quality=quals)

    expect_identical(nchar(out$quality), nchar(gsub("-", "", out$read)))
    expect_identical(nchar(out$quality), out$end - out$start + 1L)
    expect_true(all(grepl("^7.*7$", out$quality)))
})
