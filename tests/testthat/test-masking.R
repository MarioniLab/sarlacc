# This tests various aspects of the masking machinery.
# library(sarlacc); library(testthat); source("test-masking.R")

test_that("masking of bad bases works correctly", {
    CHECKFUN <- function(incoming, qualities, threshold) {          
        masked <- .Call(sarlacc:::cxx_mask_bad_bases, incoming, 
            as.list(as(qualities, "NumericList")),
            threshold)

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
    set.seed(1000)
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
})


test_that("unmasking of previously masked bases works correctly", {
    MASKFUN <- function(x) {
        DNAStringSet((gsub("[ACTG]", "N", x)))
    }
    ORIGINAL <- function(x) {
        DNAStringSet(gsub("-", "", x))
    }
    unmask_bases <- function(masked, original) {
        .Call(sarlacc:::cxx_unmask_bases, masked, original)
    }

    # Simple case, no deletions.
    seq.in <- c("acacgtagtgtcagtctaacatcagctacgttacat", # no mask
                "ACACGTAGTGTCAGTCTAACATCAGCTACGTTACAT", # all masked 
                "acacgtcgtgtcagtctaaCatcagctacgttacat", # mask in the middle
                "Acacgttgtgtcagtctaacatcagctacgttacat", # mask at the start 
                "acacgtggtgtcagtctaacatcagctacgttacaT") # mask at the end
    expect_equal(unmask_bases(MASKFUN(seq.in), ORIGINAL(seq.in)),
                 toupper(seq.in))

    # Deletions in the middle.
    seq.in <- c("acacgtagtgtcagtc-taacatcagctacgttacat",
                "ACACGTAGTGTCAGTC-TAACATCAGCTACGTTACAT",
                "acacgtcgtgtcagtc-taaCatcagctacgttacat", 
                "Acacgttgtgtcagtc-taacatcagctacgttacat", 
                "acacgtggtgtcagtc-taacatcagctacgttacaT") 
    expect_equal(unmask_bases(MASKFUN(seq.in), ORIGINAL(seq.in)),
                 toupper(seq.in))

    # Deletions at the start.
    seq.in <- c("-acacgtcgtgtcagtctaacatcagctacgttacat",
                "-ACACGTAGTGTCAGTCTAACATCAGCTACGTTACAT",
                "-acacgtcgtgtcagtctaaCatcagctacgttacat", 
                "-Acacgtcgtgtcagtctaacatcagctacgttacat", 
                "-acacgtcgtgtcagtctaacatcagctacgttacaT") 
    expect_equal(unmask_bases(MASKFUN(seq.in), ORIGINAL(seq.in)),
                 toupper(seq.in))

    # Deletions at the end.
    seq.in <- c("acacgtcgtgtcagtctaacatcagctacgttacat-",
                "ACACGTAGTGTCAGTCTAACATCAGCTACGTTACAT-",
                "acacgtcgtgtcagtctaaCatcagctacgttacat-", 
                "Acacgtcgtgtcagtctaacatcagctacgttacat-", 
                "acacgtcgtgtcagtctaacatcagctacgttacaT-") 
    expect_equal(unmask_bases(MASKFUN(seq.in), ORIGINAL(seq.in)),
                 toupper(seq.in))

    # This should trigger an error:
    expect_error(unmask_bases(DNAStringSet("AA-AA"), DNAStringSet("AAA")),
                "different lengths")
    expect_error(unmask_bases(DNAStringSet("NNNNN"), DNAStringSet("AAA")),
                "sequence in alignment string is longer than the original")
    expect_error(unmask_bases(DNAStringSet(c("AAAA", "GGGG")), DNAStringSet("AAA")),
                "alignment and original sequences should have the same number of entries")

    # This should be fine:
    expect_identical(length(unmask_bases(DNAStringSet(), DNAStringSet())), 0L)

})
