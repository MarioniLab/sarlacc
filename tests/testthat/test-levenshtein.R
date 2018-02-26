# This script tests the various levenshtein-computing machinery in sarlacc.
# library(sarlacc); library(testthat); source("test-levenshtein.R")

SEQSIM <- function(n, min.length, max.length) {
    collected <- character(n)
    for (i in seq_len(n)) {
        L <- floor(runif(1, min.length, max.length+1))
        collected[i] <- paste(sample(c("A", "C", "G", "T"), L, replace=TRUE), collapse="")
    }
    return(DNAStringSet(collected))
}

test_that("custom levenshtein calculator works as expected for default data", {
    set.seed(1000)          

    # Simulate many, many random sequences.
    randoms <- SEQSIM(100, min.length=1, max.length=20)
    ref <- stringDist(randoms, method="levenshtein")
    ref <- as.integer(ref)
    out <- .Call(sarlacc:::cxx_compute_lev_masked, randoms)
    expect_identical(ref, out) 

    # Again, with a more restricted range.
    randoms <- SEQSIM(100, min.length=5, max.length=10)
    ref <- stringDist(randoms, method="levenshtein")
    ref <- as.integer(ref)
    out <- .Call(sarlacc:::cxx_compute_lev_masked, randoms)
    expect_identical(ref, out) 
})

test_that("custom levenshtein calculator works as expected for masked data", {
    ref <- "ACAGCTAGC"
    for (i in seq_len(nchar(ref))) { 
        masked <- ref
        substr(masked, i, i+1) <- "N"
        out <- .Call(sarlacc:::cxx_compute_lev_masked, DNAStringSet(c(ref, masked)))
        expect_identical(out, 1L)

        out <- .Call(sarlacc:::cxx_compute_lev_masked, DNAStringSet(c(ref, ref)))
        expect_identical(out, 0L)
        
        # N's are missing, so a distance of 1 even when strings are the same
        out <- .Call(sarlacc:::cxx_compute_lev_masked, DNAStringSet(c(masked, masked)))
        expect_identical(out, 1L)
    }
})

test_that("fast levenshtein finder works as expected for default data", {
    set.seed(1000)          

    # Simulate many, many random sequences.
    randoms <- SEQSIM(100, min.length=1, max.length=20)
    ref <- as.matrix(stringDist(randoms, method="levenshtein"))

    closest <- .Call(sarlacc:::cxx_umi_group, randoms, order(randoms)-1L, 5L)
    for (x in seq_along(closest)) {
        expect_identical(unname(which(ref[x,] <= 5L)), sort(closest[[x]]))
    }

    # Simulate many, many random sequences.
    randoms <- SEQSIM(100, min.length=5, max.length=10)
    ref <- as.matrix(stringDist(randoms, method="levenshtein"))

    closest <- .Call(sarlacc:::cxx_umi_group, randoms, order(randoms)-1L, 5L)
    for (x in seq_along(closest)) {
        expect_identical(unname(which(ref[x,] <= 5L)), sort(closest[[x]]))
    }
})

test_that("fast levenshtein finder works as expected for masked data", {
    ref <- "ACAGCTAGC"
    for (i in seq_len(nchar(ref))) { 
        masked <- ref
        substr(masked, i, i+1) <- "N"
        
        set <- DNAStringSet(c(ref, masked))
        out <- .Call(sarlacc:::cxx_umi_group, set, order(set)-1L, 1L)
        expect_identical(sort(out[[1]]), c(1L, 2L))
        expect_identical(sort(out[[2]]), c(1L, 2L))

        # N's are missing, so a distance of 1 even when strings are the same
        out <- .Call(sarlacc:::cxx_umi_group, set, order(set)-1L, 0L)
        expect_identical(sort(out[[1]]), c(1L))
        expect_identical(sort(out[[2]]), integer(0))
    }
})
