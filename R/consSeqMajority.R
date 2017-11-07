consSeqMajority <- function(msalign, id_seq_len.fil, nullborder = 0.6)

#ConsensusSequence by majority decision
#Returns List of consensus sequences for each cluster id and a list with the corresponding phredscores (which is named by the coverage of reads per consensus)

{
    consseq.c4_maj <- vector("list", length(id_seq_len.fil))
    phredscore_c4_maj <- vector("list", length(id_seq_len.fil))
    read_coverage <- vector("list", length(id_seq_len.fil))
    
    for (i in 1:length(id_seq_len.fil)){
        aln_char <- as.character(msalign[[i]])
        aln_split <- strsplit(aln_char, "")
        x <- do.call(rbind, aln_split)
        is.A <- x=="A"
        is.T <- x=="T"
        is.C <- x=="C"
        is.G <- x=="G"
        percent.A <- colSums(is.A)/nrow(x) #counts TRUE values for each column (each position in read)
        percent.T <- colSums(is.T)/nrow(x)
        percent.C <- colSums(is.C)/nrow(x)
        percent.G <- colSums(is.G)/nrow(x)
    
        bases <- c("A","T","C","G")
        maj_base <- cbind(percent.A, percent.T, percent.C, percent.G)
    
        consseq_maj <- bases[max.col(maj_base)]
        consseq_maj[((percent.A + percent.C + percent.G + percent.T) < nullborder)] <- "-"
    
        consseq.c4_maj[[i]] <- consseq_maj
    
        phredscore <- rep("-", ncol(x))
        phredscore[consseq_maj=="A"] <- log10(1-(percent.A[consseq_maj=="A"]-0.000001))*-10
        phredscore[consseq_maj=="T"] <- log10(1-(percent.T[consseq_maj=="T"]-0.000001))*-10
        phredscore[consseq_maj=="C"] <- log10(1-(percent.C[consseq_maj=="C"]-0.000001))*-10
        phredscore[consseq_maj=="G"] <- log10(1-(percent.G[consseq_maj=="G"]-0.000001))*-10
    
        phredscore_c4_maj[[i]] <- phredscore
    
        read_coverage[[i]] <- nrow(x)
        names(phredscore_c4_maj) <- read_coverage
    
    }
    return(list(read_coverage, phredscore_c4_maj, consseq.c4_maj))
}
