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

##############################################################

test_that("fast levenshtein finder works as expected for default data", {
    set.seed(1000)          

    # Simulate many random sequences.
    randoms <- SEQSIM(100, min.length=1, max.length=20)
    ref <- as.matrix(stringDist(randoms, method="levenshtein"))

    closest <- .Call(sarlacc:::cxx_umi_group, randoms, 5L)
    for (x in seq_along(closest)) {
        expect_identical(unname(which(ref[x,] <= 5L)), sort(closest[[x]]+1L))
    }

    closest <- .Call(sarlacc:::cxx_umi_group, randoms, 2L)
    for (x in seq_along(closest)) {
        expect_identical(unname(which(ref[x,] <= 2L)), sort(closest[[x]]+1L))
    }

    # Another simulation of many random sequences.
    randoms <- SEQSIM(100, min.length=5, max.length=10)
    ref <- as.matrix(stringDist(randoms, method="levenshtein"))

    closest <- .Call(sarlacc:::cxx_umi_group, randoms, 5L)
    for (x in seq_along(closest)) {
        expect_identical(unname(which(ref[x,] <= 5L)), sort(closest[[x]]+1L))
    }

    closest <- .Call(sarlacc:::cxx_umi_group, randoms, 2L)
    for (x in seq_along(closest)) {
        expect_identical(unname(which(ref[x,] <= 2L)), sort(closest[[x]]+1L))
    }

    # Simulations involving lots of the same sequence.
    randoms <- SEQSIM(50, min.length=5, max.length=10)
    randoms <- randoms[sample(50, 100, replace=TRUE),]
    ref <- as.matrix(stringDist(randoms, method="levenshtein"))

    closest <- .Call(sarlacc:::cxx_umi_group, randoms, 5L)
    for (x in seq_along(closest)) {
        expect_identical(unname(which(ref[x,] <= 5L)), sort(closest[[x]]+1L))
    }

    closest <- .Call(sarlacc:::cxx_umi_group, randoms, 2L)
    for (x in seq_along(closest)) {
        expect_identical(unname(which(ref[x,] <= 2L)), sort(closest[[x]]+1L))
    }
})

test_that("fast levenshtein finder works as expected for masked data", {
    ref <- "ACAGCTAGC"
    for (i in seq_len(nchar(ref))) { 
        masked <- ref
        substr(masked, i, i+1) <- "N"
        
        set <- DNAStringSet(c(ref, masked))
        out <- .Call(sarlacc:::cxx_umi_group, set, 1L)
        expect_identical(sort(out[[1]]+1L), c(1L, 2L))
        expect_identical(sort(out[[2]]+1L), c(1L, 2L))

        # N's are missing, so a distance of 1 even when strings are the same
        out <- .Call(sarlacc:::cxx_umi_group, set, 0L)
        expect_identical(sort(out[[1]]+1L), c(1L))
        expect_identical(sort(out[[2]]+1L), integer(0))
    }
})

##############################################################

test_that("graph clusterer works as expected for a given list", {
    set.seed(1000)          

    for (nuniverse in c(20, 50, 100)) {
        for (nlinks in c(50, 100, 200)) {
            # Simulate a list.
            from <- sample(nuniverse, nlinks, replace=TRUE)
            to <- sample(nuniverse, nlinks, replace=TRUE)
            from <- factor(from, levels=seq_len(nuniverse))
            links <- split(to - 1L, from, drop=FALSE) # expects zero-indexed.

            # Create descending clusters. 
            clusters <- .Call(sarlacc:::cxx_descending_graph_cluster, links)
            expect_true(all(tabulate(clusters) >0L))
           
            # Checking all values with some tests.
            all.sizes <- lengths(links)

            pruned.links <- links
            for (i in seq_along(links)) {
                current <- links[[i]] + 1L
                keep <- all.sizes[current] <= all.sizes[i]
                pruned.links[[i]] <- current[keep]
            }

            G <- igraph::make_graph(rbind(rep(seq_along(pruned.links), lengths(pruned.links)), unlist(pruned.links)), n=nuniverse)
            for (chosen in split(seq_along(clusters), clusters)) {
                subG <- igraph::induced_subgraph(G, chosen)
                expect_identical(igraph::count_components(subG), 1)
            }
        }
    }
})
