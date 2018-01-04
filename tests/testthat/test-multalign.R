# This tests various aspects of the multiple alignment and consensus calculation machinery.
# library(sarlacc); library(testthat); source("test-multalign.R")

test_that("masking of bad bases works correctly", {
    seq.in <- c("AAAATTTTCCCCGGGG",
                "GGGGTTTTCCCCAAAA",
                "AAAACCCCTTTTGGGG")
    qualities <- c(PhredQuality(rep(10^c(-1, -2, -3, -4), each=4)),
                   PhredQuality(rep(10^c(-4, -3, -2, -1)/2, each=4)),
                   PhredQuality(rep(10^c(-3, -1, -2, -4)*2, each=4)))

    # Checking with a variety of thresholds of 0.01.
    for (threshold in c(0.001, 0.01, 0.05, 0.1)) { 
        masked <- sarlacc:::.mask_bad_bases(QualityScaledDNAStringSet(seq.in, qualities), log10(threshold)*-10)
        for (i in seq_along(seq.in)) {
            to.mask <- as.numeric(qualities[i]) > threshold
            bases <- strsplit(seq.in[i], "")[[1]]
            bases[to.mask] <- "N"
            expect_identical(paste(bases, collapse=""), as.character(masked[i]))
        }
    }
})

