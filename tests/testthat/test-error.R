# This tests the errorFinder function.
# library(sarlacc); library(testthat); source("test-error.R")

CHECKFUN <- function(aln){
    # Checking for errors other than insertions.
    exploded_aln <- as.matrix(aln)
    counted <- apply(exploded_aln, 2, FUN=function(x) { table(factor(x, levels=c("A", "C", "G", "T", "-"))) })
    counted <- DataFrame(t(counted), check.names=FALSE)
    colnames(counted) <- c('A','C', 'G', 'T', "deletion")
    ref <- gsub("-","",alignedSubject(aln)[1])
    ref <- strsplit(ref, split="")[[1]]
    counted <- cbind(base=ref, counted)
        
    # Counting insertion errors.
    collected <- vector("list", length(ref)+1L)
    for (i in seq_along(aln)) {
        subject_vec <- strsplit(as.character(alignedSubject(aln[i])), split="")[[1]]
        not_gap <- subject_vec!="-"
        assignments <- length(ref) - rev(cumsum(rev(not_gap))) + 1L
        insertions <- tabulate(assignments[!not_gap], length(collected))
        for (j in seq_along(insertions)) {
            collected[[j]] <- c(collected[[j]], unname(insertions[j]))
        }
    }

    collected <- lapply(collected, sort)
    collected <- lapply(collected, Rle)

    counted <- rbind(counted, DataFrame(lapply(counted[1,], FUN=function(x) { NA })))
    counted$insertion <- as(collected, "RleList")
    return(counted)
}

test_that("errorFinder works correctly", {
    # Works correctly with substitution-only errors.
    reads <- c("acgactagcacgtcagta", "acgactagcacTtcagta", "Gcgactagcacgtcagta", "acgactagcacgtcagtC", "GcgacCagcGGgtcaCCC")
    refs  <- c("acgactagcacgtcagta", "acgactagcacgtcagta", "acgactagcacgtcagta", "acgactagcacgtcagta", "acgactagcacgtcagta")

    for (x in seq_along(reads)) {
        aln <- PairwiseAlignmentsSingleSubject(pattern=DNAString(reads[x]), subject=DNAString(refs[x]), type="global")
        expect_identical(CHECKFUN(aln), errorFinder(aln)$full)
    }

    # Works correctly with deletions.
    reads <- c("acgactagcac-tcagta", "acgactagcac--cagta", "----ctagcacgtcagta", "acgactagcacgtca---", "-cgac--gc--gtca---")
    refs  <- c("acgactagcacgtcagta", "acgactagcacgtcagta", "acgactagcacgtcagta", "acgactagcacgtcagta", "acgactagcacgtcagta")

    for (x in seq_along(reads)) {
        aln <- PairwiseAlignmentsSingleSubject(pattern=DNAString(reads[x]), subject=DNAString(refs[x]), type="global")
        expect_identical(CHECKFUN(aln), errorFinder(aln)$full)
    }

    # Works correctly with insertions.
    reads <- c("acgactagcacgtcagta", "acgactagcacgtcagta", "acgacctagcacgtcagta", "acgactagcacgtcacga", "acgactggcttttcaatt")
    refs  <- c("acgactagcac-tcagta", "acgactagcac--cagta", "----cctagcacgtcagta", "acgactagcacgtca---", "--gac--gca--tcag--")

    for (x in seq_along(reads)) {
        aln <- PairwiseAlignmentsSingleSubject(pattern=DNAString(reads[x]), subject=DNAString(refs[x]), type="global")
        expect_identical(CHECKFUN(aln), errorFinder(aln)$full)
    }


    # Works correctly with a whole bundle of errors.
    reads <- c("ac--cAagcacgtcaCta", "--gaGGagcacgtcaCCa", "acgacctCCcacgtGa---", "--gacCCCcacgtTTcga", "acGactggcCtttc-att")
    refs  <- c("acgactagcac-tcagta", "acgactagcac--cagta", "----cctagcacgtcagta", "acgactagcacgtca---", "--gac--gca--tcag--")

    for (x in seq_along(reads)) {
        aln <- PairwiseAlignmentsSingleSubject(pattern=DNAString(reads[x]), subject=DNAString(refs[x]), type="global")
        expect_identical(CHECKFUN(aln), errorFinder(aln)$full)
    }

    # Works correctly in a live-fire excercise with a single subject.
    aln <- pairwiseAlignment(subject=DNAString(c("GGAAACGATCAGCTACGAACACT")), 
                             pattern=DNAStringSet(c("GGAACGTCAGCGGTACGAAACACTAAAA",
                                                    "GGAACGTCAGCGGTACGAAACACT",
                                                    "GGGGGAAACGATCAGCTACGAACACT",
                                                    "ACGATCAGCTACGAACACT")))
    out <- errorFinder(aln)
    expect_identical(CHECKFUN(aln), out$full)

    # Checking transition matrix margins.
    expect_equal(colSums(out$transition), unlist(lapply(out$full[,c("A", "C", "G", "T")], sum, na.rm=TRUE)))
    by.base <- split(out$full[,c("A", "C", "G", "T")], out$full$base)
    by.base <- lapply(by.base, FUN=function(x) { sum(as.matrix(x)) })
    expect_equal(rowSums(out$transition), unlist(by.base))
})

test_that("errorFinder fails correctly", {
    aln <- pairwiseAlignment(subject=DNAStringSet(c("GGAAACGATCAGCTACGAACACT", "GGAAACGATCAGCTACGAACACT")),
                             pattern=DNAStringSet(c("GGAACGTCAGCGGTACGAAACACTAAAA", "GGAACGTCAGCGGTACGAAACACTAAAA")))
    expect_error(errorFinder(aln), "single subject")

    aln <- pairwiseAlignment(subject=DNAStringSet(c("GGAAACGATCAGCTACGAACACT")), type="local",
                             pattern=DNAStringSet(c("GGAACGTCAGCGGTACGAAACACTAAAA", "GGAACGTCAGCGGTACGAAACACTAAAA")))
    expect_error(errorFinder(aln), "global")

    aln <- pairwiseAlignment(subject="GGAAACGATCAGCTACGAACACT", 
                             pattern="GGAACGTCAGCGGTACGAAACACTAAAA")
    expect_error(errorFinder(aln), "DNAString")
})

