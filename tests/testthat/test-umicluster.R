# This script tests the UMI-based clustering in sarlacc.
# library(sarlacc); library(testthat); source("test-umicluster.R")

REF <- function(groups)
# Slow reference function in R.
{
    collected <- vector("list", length(groups))

    for (i in seq_along(groups)){
        groupsize <- lengths(groups)
        chosen <- max(which(groupsize==max(groupsize))) # last, if ties.
        curgroup <- groups[[chosen]]
        collected[[i]] <- curgroup

        # Removing all elements in the current group from all other groups.
        for (j in curgroup) {
            groups[[j]] <- character(0)
        }
        for (j in seq_along(groups)) {
            groups[[j]] <- setdiff(groups[[j]], curgroup)
        }

        if (all(length(groups)==0)) {
            break
        }
    }

    return(collected[lengths(collected)>0])
}

library(Matrix)
MOCKUP <- function(nnodes, density)
# Mock function to generate the desired symmetric link set.
{
    out <- rsparsematrix(nnodes, nnodes, density=density, symmetric=TRUE)
    diag(out) <- 1
    has.link <- out!=0

    idx <- which(has.link, arr.ind=TRUE)
    split(idx[,1], idx[,2])
}

COMPARE <- function(left, right)
# As ordering of entries may not be the same, due to special handling of solo groups.
{
    expect_identical(length(left), length(right))
    expect_true(all(left %in% right))
    expect_true(all(right %in% left))
}

test_that("UMI clustering routine works as expected", {
    links <- MOCKUP(20, density=0.05)
    ref <- REF(links)
    obs <- .Call(sarlacc:::cxx_cluster_umis_test, links)
    COMPARE(ref, obs)

    links <- MOCKUP(20, density=0.1)
    ref <- REF(links)
    obs <- .Call(sarlacc:::cxx_cluster_umis_test, links)
    COMPARE(ref, obs)

    links <- MOCKUP(20, density=0.2)
    ref <- REF(links)
    obs <- .Call(sarlacc:::cxx_cluster_umis_test, links)
    COMPARE(ref, obs)

    links <- MOCKUP(50, density=0.2)
    ref <- REF(links)
    obs <- .Call(sarlacc:::cxx_cluster_umis_test, links)
    COMPARE(ref, obs)

    links <- MOCKUP(50, density=0.4)
    ref <- REF(links)
    obs <- .Call(sarlacc:::cxx_cluster_umis_test, links)
    COMPARE(ref, obs)

    links <- MOCKUP(50, density=0.1)
    ref <- REF(links)
    obs <- .Call(sarlacc:::cxx_cluster_umis_test, links)
    COMPARE(ref, obs)

    # Handles solo-only groups correctly.
    links <- MOCKUP(50, density=0)
    ref <- REF(links)
    obs <- .Call(sarlacc:::cxx_cluster_umis_test, links)
    COMPARE(ref, obs)
})
