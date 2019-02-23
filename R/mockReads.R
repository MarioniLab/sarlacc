#' @importFrom Biostrings reverseComplement DNAString DNAStringSet PhredQuality writeXStringSet QualityScaledDNAStringSet quality
#' @importClassesFrom Biostrings BStringSet
#' @importFrom stats runif rbinom
#' @importFrom methods as
mockReads <- function(adaptor1, adaptor2, filepath,
    all.barcodes=NULL, barcode.position="auto", umi.position="auto", 
    nmolecules=10, nreads.range=c(10, 50), seqlen.range=c(500, 5000),
    sub.rate=0.05, indel.rate=0.01, max.insert=5, flip.strands=TRUE)
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

    if (length(starts)==1L && starts==-1L) {
        barcode.position <- umi.position <- c(1L, 0L)
    } else {
        ends <- starts + attr(allNs, "match.length") - 1L

        if (identical(barcode.position, "auto")) {
            barcode.position <- c(starts[1], ends[1])
        }

        if (length(starts)==1L) {
            umi.position <- c(1L, 0L)
        } else {
            if (identical(umi.position, "auto")) {
                umi.position <- c(starts[2], ends[2])
            }
        }
    }

    barcode.len <- barcode.position[2] - barcode.position[1] + 1L
    umi.len <- umi.position[2] - umi.position[1] + 1L

    # Checking barcode content, or synthesizing them.
    if (is.null(all.barcodes)) { 
        all.barcodes <- strrep(nucleotides, barcode.len)
    } else if (!all(nchar(all.barcodes)==barcode.len)) {
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
        
        # Introducing mutations or deletions to each read.
        readseq <- quals <- vector("list", nreads)
        for (j in seq_len(nreads)) {
            reref <- ref
    
            # Substitutions.
            chosen <- rbinom(length(ref), 1, sub.rate)==1
            reref[chosen] <- sample(nucleotides, sum(chosen), replace=TRUE)
    
            # Indels.
            chosen <- rbinom(length(ref), 1, indel.rate)==1
            new.values <- strrep(reref[chosen], sample(c(0L, 2:max.insert), sum(chosen), replace=TRUE))
            reref[chosen] <- new.values
   
            readseq[[j]] <- paste(reref, collapse="")
            quals[[j]] <- PhredQuality(runif(nchar(readseq[[j]]), 0, sub.rate + indel.rate)) # Made up qualities!
        }
    
        reads <- DNAStringSet(unlist(readseq))
        qreads <- QualityScaledDNAStringSet(reads, do.call(c, quals))
        names(qreads) <- paste0("MOLECULE_", i, ":READ_", seq_len(nreads))
    
        # Reversing 50% of them.
        if (flip.strands) {
            flip <- rbinom(nreads, 1, 0.5)==1
            qreads[flip] <- reverseComplement(qreads[flip])
        }
        writeXStringSet(qreads, filepath=filepath, append=(i!=1), format="fastq", qualities=as(quality(qreads), "BStringSet"))
    }
   
    reference <- DNAStringSet(unlist(reference))
    names(reference) <- paste0("MOLECULE_", seq_along(reference))
    reference
}
