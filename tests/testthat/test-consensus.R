# This script tests various aspects of the consensus sequence machinery.
# library(sarlacc); library(testthat); source("test-consensus.R")

test.align <- c("AAAAGAAAAA-AAATAAAA",
                "ACACA-AAAA--AAT-AGA",
                "GA-AG-C-A-T-AAT-AAA",
                "AT-AG-T-AGTAAGA-AGA",
                "-AAAGAT-AGTCAGA-AGA",
                "AGAAAAT-AGAAATA-AGA") # variety of errors; indels and/or substitutions

n.align <- c("NAAAAANNN",
             "NNAANA---",
             "NNNANNN--",
             "NNNN--NN-")
    
library(Biostrings)

#####################################################

nucleotides <- c("A", "C", "G", "T")
BASICFUN <- function(current, min.coverage, pseudo.count) {
    x <- consensusMatrix(current) 
    x <- x[c(nucleotides, "N"),,drop=FALSE] 
    base.exists <- colSums(x) >= min.coverage * length(current)
    x <- x[nucleotides,base.exists,drop=FALSE]
    
    # Defining the consensus sequence based on relative majority coverage
    # (at least one base must be above 25%, as we only keep 'x' in "ATCG")
    x <- t(x)
    majority.base <- max.col(x, ties.method="first")
    conseq <- nucleotides[majority.base]
    consensus <- paste(conseq, collapse="")
    
    # Computing the error probabilities.
    n.chosen <- numeric(length(conseq))
    for (j in nucleotides) { 
        chosen <- conseq==j
        n.chosen[chosen] <- x[chosen,j]
    }
    error <- 1 - (n.chosen + pseudo.count/4)/(rowSums(x)+pseudo.count)
    return(list(consensus, log(error)))
}

test_that("basic consensus formation works properly", {
    current <- DNAStringSet(test.align)
    expect_identical(BASICFUN(current, min.coverage=0.6, pseudo.count=1),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.6, 1))
    expect_identical(BASICFUN(current, min.coverage=0.6, pseudo.count=2),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.6, 2))
    expect_identical(BASICFUN(current, min.coverage=0.2, pseudo.count=2),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.2, 2))
    expect_identical(BASICFUN(current, min.coverage=0.9, pseudo.count=1),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.9, 1))

    # Throwing in variable numbers of Ns, with and without deletions.
    # (note that N's default to A's, but obviously at very low quality).
    current <- DNAStringSet(n.align)
    expect_identical(BASICFUN(current, min.coverage=0.6, pseudo.count=1),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.6, 1))
    expect_identical(BASICFUN(current, min.coverage=0.6, pseudo.count=2),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.6, 2))
    expect_identical(BASICFUN(current, min.coverage=0.2, pseudo.count=2),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.2, 2))
    expect_identical(BASICFUN(current, min.coverage=0.9, pseudo.count=1),
                     .Call(sarlacc:::cxx_create_consensus_basic, current, 0.9, 1))

    # Returns something with an empty input.
    expect_identical(.Call(sarlacc:::cxx_create_consensus_basic, current[0], 0.6, 1), list("", numeric(0)))
})

errorToPhred <- function(errors) {
    encoding <- encoding(PhredQuality(""))
    score <- round(errors/log(10) * -10)
    score <- pmax(min(encoding), pmin(score, max(encoding)))
    m <- match(score, encoding)
    paste(names(encoding)[m], collapse="")
}

test_that("basic looped consensus formation works properly", {
    stuff <- list(DNAStringSet(n.align), DNAStringSet(test.align))
    out <- .Call(sarlacc:::cxx_create_consensus_basic_loop, stuff, 0.6, 1)

    ref1 <- .Call(sarlacc:::cxx_create_consensus_basic, stuff[[1]], 0.6, 1)
    ref2 <- .Call(sarlacc:::cxx_create_consensus_basic, stuff[[2]], 0.6, 1)
    expect_identical(out[[1]], c(ref1[[1]], ref2[[1]]))

    expect_identical(out[[2]], c(errorToPhred(ref1[[2]]), errorToPhred(ref2[[2]])))
})

#####################################################

QUALFUN <- function(alignments, qualities, min.coverage) {
    aln.mat <- do.call(rbind, strsplit(as.character(alignments), ""))

    # Forming a quality matrix.
    qual.mat <- matrix(NA_integer_, nrow(aln.mat), ncol(aln.mat))
    for (i in seq_len(nrow(aln.mat))) { 
        current <- aln.mat[i,]
        chosen <- which(current!="-")

        # Passaging through a string to mimic loss of precision.
        quals <- as(qualities[[i]], "PhredQuality")
        quals <- as.numeric(quals) 
        qual.mat[i,chosen] <- quals
    }

    # Subsetting based on the minimum coverage.
    keep <- colSums(aln.mat!="-")/nrow(aln.mat) >= min.coverage
    aln.mat <- aln.mat[,keep,drop=FALSE]
    qual.mat <- qual.mat[,keep,drop=FALSE]

    # Running through the alignment matrix.
    consensus <- character(ncol(aln.mat))
    out.err <- numeric(ncol(aln.mat))
    for (i in seq_len(ncol(aln.mat))) { 
        curcol <- aln.mat[,i]
        keep <- curcol %in% nucleotides
        curcol <- curcol[keep]

        curqual <- qual.mat[keep,i]
        correct <- log(1 - curqual)
        incorrect <- log(curqual/3)

        # Computing the probability of the current position being a particular base. 
        all.probs <- c(A=sum(correct[curcol=="A"]) + sum(incorrect[curcol!="A"]),
                       C=sum(correct[curcol=="C"]) + sum(incorrect[curcol!="C"]),
                       G=sum(correct[curcol=="G"]) + sum(incorrect[curcol!="G"]),
                       T=sum(correct[curcol=="T"]) + sum(incorrect[curcol!="T"]))
        all.probs <- exp(all.probs)
        all.probs <- all.probs/sum(all.probs)

        chosen <- which.max(all.probs)
        consensus[i] <- names(all.probs)[chosen]
        out.err[i] <- log(sum(all.probs[-chosen]))  # not using 1-chosen for numerical precision purposes.
    }
    
    return(list(paste(consensus, collapse=""), out.err))
}

.generate_qualities <- function(strings, upper) {
    qualities <- vector("list", length(strings))
    for (i in seq_along(qualities)) { 
        qualities[[i]] <- runif(nchar(gsub("-", "", strings[i])), 0, upper)
    }
    return(qualities)
}

.qualities_to_phred <- function(qualities) {
    if (length(qualities)) {
        quals <- do.call(c, lapply(qualities, as, "PhredQuality"))
    } else {
        quals <- PhredQuality(1)[0]
    }
    quals
}

TESTFUN <- function(seq, qualities, min.coverage) {
    quals <- .qualities_to_phred(qualities)
    enc <- sarlacc:::.create_encoding_vector(quals)
    .Call(sarlacc:::cxx_create_consensus_quality, seq, min.coverage, quals, enc)
}

set.seed(1000)
test_that("quality-based formation of a consensus sequence works correctly", {
    current <- DNAStringSet(test.align)
    for (upper in c(0.01, 0.1, 0.5, 1)) {
        qualities <- .generate_qualities(test.align, upper)
        expect_equal(QUALFUN(current, qualities, min.coverage=0.6), TESTFUN(current, qualities, min.coverage=0.6), tol=1e-5)
        expect_equal(QUALFUN(current, qualities, min.coverage=0.2), TESTFUN(current, qualities, min.coverage=0.2), tol=1e-5)
        expect_equal(QUALFUN(current, qualities, min.coverage=0.9), TESTFUN(current, qualities, min.coverage=0.9), tol=1e-5)
    }

    # Same sort of treatment with the N's.
    current <- DNAStringSet(n.align)
    for (upper in c(0.01, 0.1, 0.5, 1)) { 
        qualities <- vector("list", length(current))
        for (i in seq_along(qualities)) { 
            qualities[[i]] <- runif(nchar(gsub("-", "", n.align[i])), 0, upper)
        }
        expect_equal(QUALFUN(current, qualities, min.coverage=0.6), TESTFUN(current, qualities, min.coverage=0.6), tol=1e-5)
        expect_equal(QUALFUN(current, qualities, min.coverage=0.2), TESTFUN(current, qualities, min.coverage=0.2), tol=1e-5)
        expect_equal(QUALFUN(current, qualities, min.coverage=0.9), TESTFUN(current, qualities, min.coverage=0.9), tol=1e-5)
    }

    # Returns something with an empty input.
    expect_identical(TESTFUN(current[0], qualities[0], min.coverage=0.6), list("", numeric(0)))

    # Throws the proper errors.
    expect_error(TESTFUN(current[1], qualities[0], min.coverage=0.6), "different numbers of entries") 
    expect_error(TESTFUN(current[1], list(numeric(1)), min.coverage=0.6), "quality vector is shorter than the alignment sequence")
    expect_error(TESTFUN(current[1], list(runif(1000)), min.coverage=0.6), "quality vector is longer than the alignment sequence")
})

test_that("looped quality consensus formation works properly", {
    stuff <- list(test.align, n.align)
    all.quals <- lapply(stuff, .generate_qualities, upper=0.1)
    as.phred <- lapply(all.quals, .qualities_to_phred)
    out <- .Call(sarlacc:::cxx_create_consensus_quality_loop, stuff, 0.6, as.phred, sarlacc:::.create_encoding_vector(as.phred[[1]]))

    ref <- mapply(TESTFUN, seq=stuff, qualities=as.phred, MoreArgs=list(min.coverage=0.6), SIMPLIFY=FALSE)
    expect_identical(out[[1]], vapply(ref, "[[", i=1, FUN.VALUE="character"))
    expect_identical(out[[2]], c(errorToPhred(ref[[1]][[2]]), errorToPhred(ref[[2]][[2]])))
})

