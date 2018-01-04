# This tests various aspects of the multiple alignment and consensus calculation machinery.
# library(sarlacc); library(testthat); source("test-multalign.R")

test_that("masking of bad bases works correctly", {
    CHECKFUN <- function(incoming, qualities, threshold) {          
        masked <- sarlacc:::.mask_bad_bases(QualityScaledDNAStringSet(incoming, qualities), log10(threshold)*-10)
        for (i in seq_along(incoming)) {
            to.mask <- as.numeric(qualities[i]) > threshold
            bases <- strsplit(incoming[i], "")[[1]]
            bases[to.mask] <- "N"
            expect_identical(paste(bases, collapse=""), as.character(masked[i]))
        }
        invisible(NULL)
    }

    # Checking for a toy example.
    seq.in <- c("AAAATTTTCCCCGGGG",
                "GGGGTTTTCCCCAAAA",
                "AAAACCCCTTTTGGGG")
    qualities <- c(PhredQuality(rep(10^c(-1, -2, -3, -4), each=4)),
                   PhredQuality(rep(10^c(-4, -3, -2, -1)/2, each=4)),
                   PhredQuality(rep(10^c(-3, -1, -2, -4)*2, each=4)))

    for (threshold in c(0.001, 0.01, 0.05, 0.1)) { 
        CHECKFUN(seq.in, qualities, threshold)
    }

    # Repeating with some more random sequences.
    all.seq <- character(50)
    all.qual <- vector("list", length(all.seq))
    for (i in seq_along(all.seq)) {
        len <- round(runif(1, 10, 50))
        all.seq[i] <- paste(sample(c("A", "C", "T", "G"), len, replace=TRUE), collapse="")
        all.qual[[i]] <- PhredQuality(runif(len))
    }
    all.qual <- do.call(c, all.qual)
    for (threshold in c(0.001, 0.01, 0.05, 0.1)) { 
        CHECKFUN(all.seq, all.qual, threshold)
    }

    # Check that it throws with silly inputs.
    expect_error(sarlacc:::.mask_bad_bases(QualityScaledDNAStringSet("AAA", PhredQuality(0)), 20),
                 "sequence and quality strings are not the same length")
})


test_that("unmasking of previously masked bases works correctly", {
    MASKFUN <- function(x) {
        DNAStringSet((gsub("[ACTG]", "N", x)))
    }
    ORIGINAL <- function(x) {
        DNAStringSet(gsub("-", "", x))
    }

    # Simple case, no deletions.
    seq.in <- c("acacgtagtgtcagtctaacatcagctacgttacat",
                "acacgtcgtgtcagtctaaCatcagctacgttacat", # mask in the middle
                "Acacgttgtgtcagtctaacatcagctacgttacat", # mask at the start 
                "acacgtggtgtcagtctaacatcagctacgttacaT") # mask at the end
    expect_equal(sarlacc:::.unmask_bases(MASKFUN(seq.in), ORIGINAL(seq.in)),
                 DNAStringSet(seq.in))

    # Deletions in the middle.
    seq.in <- c("acacgtagtgtcagtc-taacatcagctacgttacat",
                "acacgtcgtgtcagtc-taaCatcagctacgttacat", 
                "Acacgttgtgtcagtc-taacatcagctacgttacat", 
                "acacgtggtgtcagtc-taacatcagctacgttacaT") 
    expect_equal(sarlacc:::.unmask_bases(MASKFUN(seq.in), ORIGINAL(seq.in)),
                 DNAStringSet(seq.in))

    # Deletions at the start.
    seq.in <- c("-acacgtcgtgtcagtctaacatcagctacgttacat",
                "-acacgtcgtgtcagtctaaCatcagctacgttacat", 
                "-Acacgtcgtgtcagtctaacatcagctacgttacat", 
                "-acacgtcgtgtcagtctaacatcagctacgttacaT") 
    expect_equal(sarlacc:::.unmask_bases(MASKFUN(seq.in), ORIGINAL(seq.in)),
                 DNAStringSet(seq.in))

    # Deletions at the end.
    seq.in <- c("acacgtcgtgtcagtctaacatcagctacgttacat-",
                "acacgtcgtgtcagtctaaCatcagctacgttacat-", 
                "Acacgtcgtgtcagtctaacatcagctacgttacat-", 
                "acacgtcgtgtcagtctaacatcagctacgttacaT-") 
    expect_equal(sarlacc:::.unmask_bases(MASKFUN(seq.in), ORIGINAL(seq.in)),
                 DNAStringSet(seq.in))

    # This should trigger an error:
    expect_error(sarlacc:::.unmask_bases(DNAStringSet("AA-AA"), DNAStringSet("AAA")),
                "different lengths")
    expect_error(sarlacc:::.unmask_bases(DNAStringSet("NNNNN"), DNAStringSet("AAA")),
                "sequence in alignment string is longer than the original")
    expect_error(sarlacc:::.unmask_bases(DNAStringSet(c("AAAA", "GGGG")), DNAStringSet("AAA")),
                "alignment and original sequences should have the same length")

    # This should be fine:
    expect_identical(length(sarlacc:::.unmask_bases(DNAStringSet(), DNAStringSet())), 0L)

})
