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

    idx <- Matrix::which(has.link, arr.ind=TRUE)
    split(idx[,1], idx[,2])
}

COMPARE <- function(left, right)
# As ordering of entries may not be the same, due to special handling of solo groups.
{
    expect_identical(length(left), length(right))
    expect_true(all(left %in% right))
    expect_true(all(right %in% left))
}

####################################################################

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

####################################################################

library(testthat); library(sarlacc)
SEQSIM <- function(n, len) {
    collected <- character(n)
    ref <- sample(c("A", "C", "G", "T"), len, replace=TRUE)
    for (i in seq_len(n)) {
        tmp <- ref
        chosen <- rbinom(len, 1, 0.1)==1 
        tmp[chosen] <- sample(c("A", "C", "G", "T"), sum(chosen), replace=TRUE)
        collected[i] <- paste(tmp, collapse="")
    }
    return(DNAStringSet(collected))
}

REFVERSION <- function(UMI1, threshold1 = 3, UMI2 = NULL, threshold2 = threshold1) { 
    out1 <- .Call(sarlacc:::cxx_fast_levdist_test, UMI1, threshold1, TRUE)

    # Repeating for the second UMI, if it is available.
    if (!is.null(UMI2)) {
        if (length(UMI1)!=length(UMI2)) { 
            stop("mismatch in lengths between 'UMI1' and 'UMI2'")
        }
        out2 <- .Call(sarlacc:::cxx_fast_levdist_test, UMI2, threshold2, TRUE)
        out1 <- mapply(intersect, out2, out1, SIMPLIFY=FALSE) # out2 order, to match C++ umi group code.
    }

    .Call(sarlacc:::cxx_cluster_umis_test, out1)
}

test_that("UMI grouping works with one UMI", {
    seqs <- list()
    for (x in 1:10) {
        N <- round(runif(1, 10, 20))
        seqs[[x]] <- SEQSIM(N, 10)
    }

    pregroups <- rep(seq_along(seqs), lengths(seqs))
    seqs <- do.call(c, seqs)

    # Shuffling to provide a more robust test.
    o <- sample(length(pregroups))
    pregroups <- pregroups[o]
    seqs <- DNAStringSet(seqs[o])

    # Testing umiGroup against a reference implementation.
    obs <- umiGroup(seqs, threshold1=1)
    ref <- REFVERSION(seqs, threshold1=1)
    expect_equal(obs, ref)

    obs <- umiGroup(seqs, threshold1=3)
    ref <- REFVERSION(seqs, threshold1=3)
    expect_equal(obs, ref)

    # With groups.
    obs <- umiGroup(seqs, threshold1=1, groups=pregroups)
    by.group <- split(seq_along(pregroups), pregroups)
    alt <- umiGroup(seqs, threshold1=1, groups=by.group)
    expect_equal(alt, obs)

    ref <- lapply(split(seqs, pregroups), REFVERSION, threshold1=1)
    for (x in seq_along(ref)) {
        i <- by.group[[x]]
        ref[[x]] <- lapply(ref[[x]], function(v) { i[v] })
    }
    expect_equal(obs, unlist(ref, recursive=FALSE, use.names=FALSE))
})

test_that("UMI grouping works with two UMIs", {
    seqs1 <- seqs2 <- list()
    for (x in 1:10) {
        N <- round(runif(1, 10, 20))
        seqs1[[x]] <- SEQSIM(N, 10)
        seqs2[[x]] <- SEQSIM(N, 5)
    }

    pregroups <- rep(seq_along(seqs1), lengths(seqs1))
    seqs1 <- do.call(c, seqs1)
    seqs2 <- do.call(c, seqs2)

    # Shuffling to provide a more robust test.
    o <- sample(length(pregroups))
    pregroups <- pregroups[o]
    seqs1 <- DNAStringSet(seqs1[o])
    seqs2 <- DNAStringSet(seqs2[o])

    # Testing umiGroup against a reference implementation.
    obs <- umiGroup(seqs1, threshold1=1, UMI2=seqs2, threshold2=1)
    ref <- REFVERSION(seqs1, threshold1=1, UMI2=seqs2, threshold2=1)
    expect_equal(obs, ref)

    obs <- umiGroup(seqs1, threshold1=3, UMI2=seqs2, threshold2=1)
    ref <- REFVERSION(seqs1, threshold1=3, UMI2=seqs2, threshold2=1)
    expect_equal(obs, ref)

    # With groups.
    obs <- umiGroup(seqs1, threshold1=1, UMI2=seqs2, threshold2=1, groups=pregroups)
    by.group <- split(seq_along(pregroups), pregroups)
    alt <- umiGroup(seqs1, threshold1=1, UMI2=seqs2, threshold2=1, groups=by.group)
    expect_equal(alt, obs)

    by.umi1 <- split(seqs1, pregroups)
    by.umi2 <- split(seqs2, pregroups)
    ref <- mapply(REFVERSION, UMI1=by.umi1, UMI2=by.umi2, threshold1=1, threshold2=1)
    for (x in seq_along(ref)) {
        i <- by.group[[x]]
        ref[[x]] <- lapply(ref[[x]], function(v) { i[v] })
    }
    expect_equal(obs, unlist(ref, recursive=FALSE, use.names=FALSE))
})


