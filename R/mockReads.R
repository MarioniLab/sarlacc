#' @importFrom Biostrings reverseComplement DNAString DNAStringSet PhredQuality writeXStringSet QualityScaledDNAStringSet quality
#' @importClaasesFrom Biostrings BStringSet
#' @importFrom stats runif rbinom
#' @importFrom methods as
mockReads <- function(adaptor1, adaptor2, filepath,
        all.barcodes=NULL, barcode.position="auto", umi.position="auto", 
        nmolecules=10, nreads.range=c(10, 50), seqlen.range=c(500, 5000)) 
# Simulates reads for examples and vignettes.
# Assumes that the first adaptor has the barcode (first Ns) and UMI (second Ns).
#
# written by Aaron Lun
# created 30 September 2018
{
    rc.adaptor2 <- reverseComplement(DNAString(adaptor2))
    nucleotides <- c("A", "C", "G", "T")

    # Checking barcode and UMI position, if not supplied.
    allNs <- gregexpr("N+", adaptor1)[[1]]
    starts <- as.integer(allNs)
    ends <- starts + attr(allNs, "match.length")

    if (identical(barcode.position, "auto")) {
        barcode.position <- c(starts[1], ends[1])
    }
    barcode.len <- barcode.position[2] - barcode.position[1] + 1L

    if (identical(umi.position, "auto")) {
        umi.position <- c(starts[2], ends[2])
    }
    umi.len <- umi.position[2] - umi.position[1] + 1L

    # Checking barcode content, or synthesizing them.
    if (is.null(all.barcodes)) { 
        all.barcodes <- strrep(nucleotides, barcode.len)
    } 
    if (!all(nchar(all.barcodes)==barcode.len)) {
        stop("'barcodes' width must correspond to barcode position")
    }

    # Looping over all molecules.
    reference <- vector("list", nmolecules)
    for (i in seq_len(nmolecules)) {
        nreads <- runif(1, nreads.range[1], nreads.range[2])
        seqlen <- runif(1, seqlen.range[1], seqlen.range[2])
        refseq <- sample(nucleotides, seqlen, replace=TRUE)
    
        # Choosing a barcode and UMI.
        tmp.adaptor1 <- adaptor1
        barcode <- sample(all.barcodes, 1)
        substr(tmp.adaptor1, barcode.position[1], barcode.position[2]) <- barcode
        umi <- paste(sample(nucleotides, umi.len, replace=TRUE), collapse="")
        substr(tmp.adaptor1, umi.position[1], umi.position[2]) <- umi
    
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
    
        reads <- DNAStringSet(unlist(readseq))
        qreads <- QualityScaledDNAStringSet(reads, do.call(c, quals))
        names(qreads) <- paste0("MOLECULE_", i, ":READ_", seq_len(nreads))
    
        # Reversing 50% of them.
        flip <- rbinom(nreads, 1, 0.5)==1
        qreads[flip] <- reverseComplement(qreads[flip])
        writeXStringSet(qreads, filepath=filepath, append=(i!=1), format="fastq", qualities=as(quality(qreads), "BStringSet"))
    }
   
    reference <- DNAStringSet(unlist(reference))
    names(reference) <- paste0("MOLECULE_", seq_along(reference))
    reference
}
