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

test_that("position shift calculations are correct", {
    position <- c(11, 14) 

    # Insert deletions throughout to mimic an alignment string for the adaptor.
    additional <- c("ACAGTGACGTNNNNACACT",
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
    test <-  sarlacc:::.compute_position_bump(additional, position)
    expect_identical(SHIFTFUN(additional, position), test)

    # Checking that they're both actually correct.    
    umi <- substr(additional, position[1] + test$start, position[2] + test$end)
    expect_identical(gsub("-", "", umi), rep("NNNN", length(additional)))
}) 
