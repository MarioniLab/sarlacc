# This script tests various aspects of the UMI extraction machinery.
# library(sarlacc); library(testthat); source("test-extract.R")

SHIFTFUN <- function(align.str, position) {
    all.dels <- gregexpr("-", align.str)
    bump.start <- bump.end <- integer(length(align.str))
    for (i in seq_along(bump.start)) { 
        dels <- as.integer(all.dels[[i]])
        if (length(dels)==1L && dels==-1L) { 
            next
        }
        true.pos <- dels - seq_along(dels)
        bump.start[i] <- sum(true.pos < position[1])
        bump.end[i] <- sum(true.pos < position[2])
    }
    return(list(start=bump.start, end=bump.end))
}

# Insert deletions throughout to mimic an alignment string for the adaptor.
position <- c(11, 14) 
adaptor_aln <- c("ACAGTGACGTNNNNACACT",
                "ACAGTG-ACGTNNNNACACT", # before the UMI
                "ACAGTGACGTNNNNACA-CT", # after the UMI
                "ACAGTGACGTNN-NNACACT", # within the UMI... and so on.
                "ACAGTGACGTN-NN-NACACT", 
                "ACAG-TGACGTNN-NNACACT",
                "A-CAGTGAC-GTNN-NNACACT",
                "-ACAGTGACGTNNNNACACT",
                "ACAGTGACGT-NNNNACACT",
                "ACAGTGACGTNNNN-ACACT",
                "AC-AGTG-ACGT-NN-NN-ACACT")

test_that("base position shift calculations are correct", {
    test <- .Call(sarlacc:::cxx_adjust_basepos_for_gaps, adaptor_aln, position[1], position[2])
    ref <- SHIFTFUN(adaptor_aln, position)
    expect_identical(ref$start, test[[1]])
    expect_identical(ref$end, test[[2]])

    # Checking that they're both actually correct.    
    umi <- substr(adaptor_aln, position[1] + test[[1]], position[2] + test[[2]])
    expect_identical(gsub("-", "", umi), rep("NNNN", length(adaptor_aln)))
})

test_that("base position shift errors are thrown", {
    expect_error(.Call(sarlacc:::cxx_adjust_basepos_for_gaps, adaptor_aln, -1L, 0L), "positive integers")
    expect_error(.Call(sarlacc:::cxx_adjust_basepos_for_gaps, adaptor_aln, 5L, 0L), "must be less than or equal to")
    expect_error(.Call(sarlacc:::cxx_adjust_basepos_for_gaps, adaptor_aln, 1L, 100L), "exceeds the number of bases")
})

test_that("alignment position shift calculations are correct", {
    read_aln <- gsub("N", "G", adaptor_aln)
    align.pos <- SHIFTFUN(read_aln, position)
    align.start <- align.pos$start + position[1]
    align.end <- align.pos$end + position[2]

    test <- .Call(sarlacc:::cxx_adjust_alignpos_for_gaps, read_aln, align.start, align.end)
    bump.start <- bump.end <- integer(length(read_aln))
    for (i in seq_along(bump.start)) { 
        is.gap <- strsplit(read_aln[i], "")[[1]]=="-"
        bump.start[i] <- sum(is.gap[seq_len(align.start[i])])
        bump.end[i] <- sum(is.gap[seq_len(align.end[i])])
    }
    expect_identical(test[[1]], bump.start)
    expect_identical(test[[2]], bump.end)

    # Checking that they're both actually correct.    
    umi <- substr(gsub("-", "", read_aln), align.start - test[[1]], align.end - test[[2]])
    expect_identical(gsub("-", "", umi), rep("GGGG", length(read_aln)))

    # Testing what happens when the start and end lie exactly on a gap
    # (this should not, in general, happen, as adjusted base positions should lie on bases).
    test <- .Call(sarlacc:::cxx_adjust_alignpos_for_gaps, "-AA-", 1, 4)
    expect_identical(test[[1]], 0L)
    expect_identical(test[[2]], 2L)
}) 

test_that("alignment position shift errors are thrown", {
    expect_error(.Call(sarlacc:::cxx_adjust_alignpos_for_gaps, adaptor_aln, -1L, 0L), "positive integers")
    expect_error(.Call(sarlacc:::cxx_adjust_alignpos_for_gaps, adaptor_aln, 5L, 0L), "must be less than or equal to")
    expect_error(.Call(sarlacc:::cxx_adjust_alignpos_for_gaps, adaptor_aln, 1L, 100L), "exceeds the number of bases")
})


