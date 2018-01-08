set.seed(1000)
adaptor1 <- "AAAACCCCNNNNGGGGTTTT" # UMI-containing adaptor.
adaptor2 <- "AAGGCCTTTTCCGGAA"     # adaptor for PCR.

listed1 <- strsplit(adaptor1, "")[[1]]
rc.adaptor2 <- reverseComplement(DNAString(adaptor2))
listed2 <- strsplit(as.character(rc.adaptor2), "")[[1]]

nucleotides <- c("A", "C", "G", "T")

fastq <- tempfile(fileext=".fastq")
for (i in seq_len(10)) {
    nreads <- runif(1, 10, 50)
    seqlen <- runif(1, 500, 5000)
    refseq <- sample(nucleotides, seqlen, replace=TRUE)
    ref <- c(listed1, refseq, listed2)
    
    # Introducing mutations or deletions.
    readseq <- quals <- vector("list", nreads)
    for (j in seq_len(nreads)) {
        reref <- ref

        # Adding substitutions at a ~10% rate.
        chosen <- rbinom(length(ref), 1, 0.1)==1
        reref[chosen] <- sample(nucleotides, sum(chosen), replace=TRUE)

        # Adding indels at a 5% rate.
        chosen <- rbinom(length(ref), 1, 0.05)==1
        new.values <- paste0(ref[chosen], sample(nucleotides, sum(chosen), replace=TRUE))
        reref[chosen] <- new.values

        readseq[[j]] <- paste(reref, collapse="")
        quals[[j]] <- PhredQuality(runif(nchar(readseq[[j]]), 0, 0.1))
    }

    collected <- DNAStringSet(unlist(readseq))
    quals <- as(do.call(c, quals), "BStringSet")
    names(collected) <- paste0("MOLECULE_", i, ":READ_", seq_len(nreads))
    writeXStringSet(collected, filepath=fastq, append=(i!=1), 
                    format="fastq", qualities=quals)
}

