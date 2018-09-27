set.seed(1000)
adaptor1 <- "ACGCAGATCGATCGATNNNNNNNNNNNNCGCGCGAGCTGACTNNNNGCACGACTCTGGTTTTTTTTTTTT" # UMI-containing adaptor.
adaptor2 <- "AAGGCCTTTTCCGACTCATGAA"     # adaptor for PCR.
rc.adaptor2 <- reverseComplement(DNAString(adaptor2))

nucleotides <- c("A", "C", "G", "T")
all.barcodes <- c("AAACCCGGGTTT", "TTTCCCGGGAAA", "TTTGGGCCCAAA", "TTTAAACCCGGG")

fastq <- tempfile(fileext=".fastq")
reference <- list()
for (i in seq_len(10)) {
    nreads <- runif(1, 10, 50)
    seqlen <- runif(1, 500, 5000)
    refseq <- sample(nucleotides, seqlen, replace=TRUE)

    # Choosing a barcode (replace the first stretch of N's).
    tmp.adaptor1 <- adaptor1
    tmp.adaptor1 <- sub("NNNNNNNNNNNN", sample(all.barcodes, 1), tmp.adaptor1)

    # Choosing a UMI (replace the second stretch of N's).
    umi <- paste(sample(nucleotides, 4, replace=TRUE), collapse="")
    tmp.adaptor1 <- sub("NNNN", umi, tmp.adaptor1)

    ref <- c(strsplit(tmp.adaptor1, "")[[1]], refseq, strsplit(as.character(rc.adaptor2), "")[[1]])
    reference[[i]] <- paste(ref, collapse="")
    
    # Introducing mutations or deletions.
    readseq <- quals <- vector("list", nreads)
    for (j in seq_len(nreads)) {
        reref <- ref

        # Adding substitutions at a ~10% rate.
        chosen <- rbinom(length(ref), 1, 0.1)==1
        reref[chosen] <- sample(nucleotides, sum(chosen), replace=TRUE)

        # Adding indels at a 5% rate.
        chosen <- rbinom(length(ref), 1, 0.05)==1
        new.values <- strrep(reref[chosen], sample(0:3, sum(chosen), replace=TRUE))
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

reference <- DNAStringSet(unlist(reference))
names(reference) <- paste0("MOLECULE_", seq_along(reference))
