# Tests the homopolymer finding and matching machinery.
# library(sarlacc); library(testthat); source("test-homopolymer.R")

FINDCHECK <- function(seq)
# Checks homopolymerFinder using RLEs.
{
    seq <- as.character(seq)
    output <- vector("list", length(seq))    
    for (i in seq_along(seq)) {
        cur.seq <- strsplit(seq[i], "")[[1]]
        cur.seq <- cur.seq[cur.seq!="-"]
        cur.rle <- rle(cur.seq)

        ends <- cumsum(cur.rle$lengths)
        starts <- ends - cur.rle$lengths + 1L
        cur.ranges <- IRanges(start=starts, ends)
        elementMetadata(cur.ranges)$base <- cur.rle$values

        output[[i]] <- cur.ranges[width(cur.ranges)>1L]
    }

    # Somewhat contrived way to generate the exact same IRangesList.
    ilist <- split(do.call(c, output), rep(seq_along(output), lengths(output)))
    names(ilist) <- names(seq)
    return(ilist)
}

test_that("homopolymer regions are found correctly", {
    # Insertions should not affect the result
    test_seq <- DNAStringSet(c("ATGG--GGCCGTTTAA",
                               "ATG-GGGCCGTTTAA",
                               "ATGGGGCCGT-TTAA",
                               "ATGGGGCCGTTTAA-A-A-",
                               "ATGGGGC-CGTTTAA",
                               "---ATGGGGCCGTTTAA"))

    ref <- FINDCHECK(test_seq)
    out <- homopolymerFinder(test_seq)
    expect_identical(ref, out)

    names(test_seq) <- paste0("SEQUENCE", seq_along(test_seq))
    ref <- FINDCHECK(test_seq)
    out <- homopolymerFinder(test_seq)
    expect_identical(ref, out)

    # Works with empty input. 
    out <- homopolymerFinder(test_seq[0])
    expect_s4_class(out, "IRangesList")
    expect_identical(length(out), 0L)
    expect_identical(mcols(unlist(out))$base, character(0))
})

##############################################################################

.get_runs_with_gaps <- function(seq, min=1L)
# Helper function to define runs within a gapped sequence.
# Returns 'index', the position of each non-gap character in the original sequence;
# and 'runs', the RLE start/end/value information for the de-gapped sequence.
{
    non.gap <- which(seq!='-')
    non.gap.rle <- rle(seq[non.gap])
    non.gap.end0 <- cumsum(non.gap.rle$lengths)
    non.gap.start0 <- non.gap.end0 - non.gap.rle$lengths + 1L
    keep <- non.gap.rle$lengths > min
    list(index=non.gap, runs=data.frame(base=non.gap.rle$value, start=non.gap.start0, end=non.gap.end0, stringsAsFactors=FALSE)[keep,])
}

MATCHCHECK <- function(alignments) 
# Checks for homopolymer matches betweeen reference and query sequences.
{
    reads <- alignedPattern(alignments)
    ref <- alignedSubject(alignments)
    output <- homopolymerFinder(ref[1])[[1]]

    total.collected <- vector("list", length(reads))
    for (i in seq_along(reads)) {
        cur.ref <- strsplit(as.character(ref[i]), "")[[1]]
        cur.read <- strsplit(as.character(reads[i]), "")[[1]]

        # Identify the reference homopolymer runs.
        ref.rle <- .get_runs_with_gaps(cur.ref, min=1L)
        ref.indices <- ref.rle$index
        ref.runs <- ref.rle$runs

        ngaps <- nrow(ref.runs)
        collected <- integer(length(output))
        expect_identical(ngaps, length(collected))

        for (run in seq_len(ngaps)) {
            s0 <- ref.runs$start[run]
            e0 <- ref.runs$end[run]
            
            # Identify the "true" run start and end on the gapped sequence.
            # Also extend the run to include adjacent gap characters.
            run.start <- ref.indices[s0]
            run.end <- ref.indices[e0]
            pre.start <- if(s0!=1L) ref.indices[s0-1]+1L else 1L
            post.end <- if(e0!=length(ref.indices)) ref.indices[e0+1]-1L else length(cur.ref) 

            # Check the corresponding read subsequence for runs of the same base.
            sub.rle <- .get_runs_with_gaps(cur.read[pre.start:post.end], min=0L)
            sub.indices <- sub.rle$index
            sub.runs <- sub.rle$runs
            sub.runs <- sub.runs[sub.runs$base==ref.runs$base[run],]

            # Figure out which read runs overlap with the "true" reference run, and pick the longest.
            keep <- overlapsAny(IRanges(sub.indices[sub.runs$start] + pre.start - 1L, 
                                        sub.indices[sub.runs$end] + pre.start - 1L), 
                                IRanges(run.start, run.end))
            run.lengths <- sub.runs$end[keep] - sub.runs$start[keep] + 1L
            collected[run] <- if(length(run.lengths)) max(run.lengths) else 0L
        }

        total.collected[[i]] <- collected
    }

    # Finally, aggregating into a list of integer vectors.
    full.list <- vector("list", length(output))
    aggregated.stats <- do.call(rbind, total.collected)
    expect_identical(ncol(aggregated.stats), length(full.list))
    for (i in seq_len(ncol(aggregated.stats))) {
        full.list[[i]] <- Rle(sort(aggregated.stats[,i]))
    }
    mcols(output)$observed <- as(full.list, "RleList")
    return(output)
}

test_that("homopolymer regions are matched correctly", {
    # Handling indels.
    reads <- c("acgtAAAAAtgca", "acgtAA-AAAtgca", "acgtAAAAAtgca", "acgtAAAAAtgca", "acgtAAAAAtgca", "acgt-AAAAAtgca", "acgtAAAAA-tgca")
    refs  <- c("acgtAAAAAtgca", "acgtAAAAAAtgca", "acgtA-AAAtgca", "acgt-AAAAtgca", "acgtAAAA-tgca", "acgtAAAAAAtgca", "acgtAAAAAAtgca")

    for (x in seq_along(reads)) {
        aln <- PairwiseAlignmentsSingleSubject(pattern=DNAString(reads[x]), subject=DNAString(refs[x]), type="global")
        expect_identical(MATCHCHECK(aln), homopolymerMatcher(aln))
    }

    # Handling substitutions within the run. 
    reads <- c("acgtAATAAtgca", "acgtTAAAAAtgca", "acgtAAACAtgca", "acgtAAAAGtgca")
    refs  <- c("acgtAAAAAtgca", "acgtAAAAAAtgca", "acgtAAAAAtgca", "acgtAAAAAtgca")

    for (x in seq_along(reads)) {
        aln <- PairwiseAlignmentsSingleSubject(pattern=DNAString(reads[x]), subject=DNAString(refs[x]), type="global")
        expect_identical(MATCHCHECK(aln), homopolymerMatcher(aln))
    }

    # Homopolymers at the start or end, with multiple deletions.
    reads <- c("CCCCCtgca", "--CCCCtgca", "CCCCCCtgca", "CCCC--tgca", "tgcaCCCCC", "tgcaCCCC--", "tgcaCCCCCC", "tgca--CCCC")
    refs  <- c("CCCCCtgca", "CCCCCCtgca", "--CCCCtgca", "CCCCCCtgca", "tgcaCCCCC", "tgcaCCCCCC", "tgcaCCCC--", "tgcaCCCCCC")

    for (x in seq_along(reads)) {
        aln <- PairwiseAlignmentsSingleSubject(pattern=DNAString(reads[x]), subject=DNAString(refs[x]), type="global")
        expect_identical(MATCHCHECK(aln), homopolymerMatcher(aln))
    }

    # Handling multiple homopolymers within the sequence, with different error types.
    reads <- c("CCCCaGGGaTT", "CCCCaCCCaTT", "CCCCaCCCaCC", "CCCCaGG-aT-", "CCGCaGAGaTC", "---CaGGGaTT", "CCCCaGGGaTTTTT")
    refs  <- c("CCCCaGGGaTT", "CCCCaCCCaTT", "CCCCaCCCaCC", "CCCCaGGGaTT", "CCCCaGGGaTT", "CCCCaGGGaTT", "CCCCaGGGa--TT-")
        
    for (x in seq_along(reads)) {
        aln <- PairwiseAlignmentsSingleSubject(pattern=DNAString(reads[x]), subject=DNAString(refs[x]), type="global")
        expect_identical(MATCHCHECK(aln), homopolymerMatcher(aln))
    }

    # Homopolymers that are observed in the reads as one or zero bases.
    reads <- c("actg--G--tt", "actg-----tt", "----tgca", "--A-tgca", "tgca----", "tgca-T--")
    refs  <- c("actgGGGGGtt", "actgGGGGGtt", "AAAAtgca", "AAAAtgca", "tgcaTTTT", "tgcaTTTT")

    for (x in seq_along(reads)) {
        aln <- PairwiseAlignmentsSingleSubject(pattern=DNAString(reads[x]), subject=DNAString(refs[x]), type="global")
        expect_identical(MATCHCHECK(aln), homopolymerMatcher(aln))
    }

    # Homopolymers that have a different base as the majority.
    reads <- c("actgAAACtttt", "actgAA-Ctttt", "actgA-A-tttt", "actg--A-tttt")
    refs  <- c("actgCCCCtttt", "actgCCCCtttt", "actgCCCCtttt", "actgCCCCtttt")

    for (x in seq_along(reads)) {
        aln <- PairwiseAlignmentsSingleSubject(pattern=DNAString(reads[x]), subject=DNAString(refs[x]), type="global")
        expect_identical(MATCHCHECK(aln), homopolymerMatcher(aln))
    }

    # Handling multiple alignments to a single subject (using multiple gap penalties to give different alignments).
    refs  <-   "AAAACCCCCGGGGGTTTT"
    reads <- c("AAAACCCCCGGGGGTTTT",     # same
               "AAAACCCCGGGGTTT",        # deletions
               "AAAAACCCCCGGGGGGTTTTTT", # insertions
               "AATACCGCCGGTGGTTAT")     # substitutions

    aln10 <- pairwiseAlignment(pattern=DNAStringSet(reads), subject=DNAStringSet(refs), gapOpening=10, gapExtension=1, 
        substitutionMatrix=nucleotideSubstitutionMatrix(), type="global")
    expect_identical(MATCHCHECK(aln10), homopolymerMatcher(aln10))

    aln1 <- pairwiseAlignment(pattern=DNAStringSet(reads), subject=DNAStringSet(refs), gapOpening=1, gapExtension=1, 
        substitutionMatrix=nucleotideSubstitutionMatrix(), type="global")
    expect_identical(MATCHCHECK(aln1), homopolymerMatcher(aln1))

    aln0 <- pairwiseAlignment(pattern=DNAStringSet(reads), subject=DNAStringSet(refs), gapOpening=0, gapExtension=0, 
        substitutionMatrix=nucleotideSubstitutionMatrix(), type="global")
    expect_identical(MATCHCHECK(aln0), homopolymerMatcher(aln0))
})

test_that("homopolymerMatcher fails correctly", {
    aln <- pairwiseAlignment(subject=DNAStringSet(c("GGAAACGATCAGCTACGAACACT", "GGAAACGATCAGCTACGAACACT")),
                             pattern=DNAStringSet(c("GGAACGTCAGCGGTACGAAACACTAAAA", "GGAACGTCAGCGGTACGAAACACTAAAA")))
    expect_error(homopolymerMatcher(aln), "single subject")

    aln <- pairwiseAlignment(subject=DNAStringSet(c("GGAAACGATCAGCTACGAACACT")), type="local",
                             pattern=DNAStringSet(c("GGAACGTCAGCGGTACGAAACACTAAAA", "GGAACGTCAGCGGTACGAAACACTAAAA")))
    expect_error(homopolymerMatcher(aln), "global")

    aln <- pairwiseAlignment(subject="GGAAACGATCAGCTACGAACACT", 
                             pattern="GGAACGTCAGCGGTACGAAACACTAAAA")
    expect_error(homopolymerMatcher(aln), "DNAString")
})

