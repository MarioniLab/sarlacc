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

##############################################################

test_that("custom levenshtein calculator works as expected for default data", {
    set.seed(1000)          

    # Simulate many, many random sequences.
    randoms <- SEQSIM(100, min.length=1, max.length=20)
    ref <- stringDist(randoms, method="levenshtein")
    out <- .Call(sarlacc:::cxx_compute_lev_masked, randoms)
    expect_identical(as.double(ref), out) 

    # Again, with a more restricted range.
    randoms <- SEQSIM(100, min.length=5, max.length=10)
    ref <- stringDist(randoms, method="levenshtein")
    out <- .Call(sarlacc:::cxx_compute_lev_masked, randoms)
    expect_identical(as.double(ref), out) 
})

test_that("custom levenshtein calculator works as expected for masked data", {
    ref <- "ACAGCTAGC"
    for (i in seq_len(nchar(ref))) { 
        masked <- ref
        substr(masked, i, i+1) <- "N"
        out <- .Call(sarlacc:::cxx_compute_lev_masked, DNAStringSet(c(ref, masked)))
        expect_identical(out, 0.5)

        out <- .Call(sarlacc:::cxx_compute_lev_masked, DNAStringSet(c(ref, ref)))
        expect_identical(out, 0)
        
        # N's are missing, so a distance of 0.5 even when strings are the same
        out <- .Call(sarlacc:::cxx_compute_lev_masked, DNAStringSet(c(masked, masked)))
        expect_identical(out, 0.5)
    }
})

##############################################################

test_that("fast levenshtein finder works as expected for default data", {
    set.seed(1000)          

    # Simulate many random sequences.
    randoms <- SEQSIM(100, min.length=1, max.length=20)
    ref <- as.matrix(stringDist(randoms, method="levenshtein"))

    closest <- .Call(sarlacc:::cxx_fast_levdist_test, randoms, 5L, TRUE)
    for (x in seq_along(closest)) {
        expect_identical(unname(which(ref[x,] <= 5L)), sort(closest[[x]]))
    }
    expect_identical(closest, .Call(sarlacc:::cxx_fast_levdist_test, randoms, 5L, FALSE))

    closest <- .Call(sarlacc:::cxx_fast_levdist_test, randoms, 2L, TRUE)
    for (x in seq_along(closest)) {
        expect_identical(unname(which(ref[x,] <= 2L)), sort(closest[[x]]))
    }
    expect_identical(closest, .Call(sarlacc:::cxx_fast_levdist_test, randoms, 2L, FALSE))

    # Another simulation of many random sequences.
    randoms <- SEQSIM(100, min.length=5, max.length=10)
    ref <- as.matrix(stringDist(randoms, method="levenshtein"))

    closest <- .Call(sarlacc:::cxx_fast_levdist_test, randoms, 5L, TRUE)
    for (x in seq_along(closest)) {
        expect_identical(unname(which(ref[x,] <= 5L)), sort(closest[[x]]))
    }
    expect_identical(closest, .Call(sarlacc:::cxx_fast_levdist_test, randoms, 5L, FALSE))

    closest <- .Call(sarlacc:::cxx_fast_levdist_test, randoms, 2L, TRUE)
    for (x in seq_along(closest)) {
        expect_identical(unname(which(ref[x,] <= 2L)), sort(closest[[x]]))
    }
    expect_identical(closest, .Call(sarlacc:::cxx_fast_levdist_test, randoms, 2L, FALSE))
})

test_that("fast levenshtein finder works as expected for highly dense space", {
    # Simulations involving lots of the same sequence.
    randoms <- SEQSIM(50, min.length=5, max.length=10)
    randoms <- randoms[sample(50, 100, replace=TRUE),]
    ref <- as.matrix(stringDist(randoms, method="levenshtein"))

    closest <- .Call(sarlacc:::cxx_fast_levdist_test, randoms, 5L, TRUE)
    for (x in seq_along(closest)) {
        expect_identical(unname(which(ref[x,] <= 5L)), sort(closest[[x]]))
    }
    expect_identical(closest, .Call(sarlacc:::cxx_fast_levdist_test, randoms, 5L, FALSE))

    closest <- .Call(sarlacc:::cxx_fast_levdist_test, randoms, 2L, TRUE)
    for (x in seq_along(closest)) {
        expect_identical(unname(which(ref[x,] <= 2L)), sort(closest[[x]]))
    }
    expect_identical(closest, .Call(sarlacc:::cxx_fast_levdist_test, randoms, 2L, FALSE))
})

test_that("fast levenshtein finder works as expected for masked data", {
    ref <- "ACAGCTAGC"
    for (i in seq_len(nchar(ref))) { 
        masked <- ref
        substr(masked, i, i+1) <- "N"
        
        set <- DNAStringSet(c(ref, masked))
        out <- .Call(sarlacc:::cxx_fast_levdist_test, set, 1L, TRUE)
        expect_identical(sort(out[[1]]), c(1L, 2L))
        expect_identical(sort(out[[2]]), c(1L, 2L))

        # N's are missing, so a distance of 1 even when strings are the same
        out <- .Call(sarlacc:::cxx_fast_levdist_test, set, 0L, TRUE)
        expect_identical(sort(out[[1]]), c(1L))
        expect_identical(sort(out[[2]]), integer(0))
    }
})
