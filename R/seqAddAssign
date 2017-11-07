seqAddAssign <- function(fasta_file, subgroups, id_reads)

#Returns list of cluster id groups that contain the corresponding sequences

{
    names.cell4 <- read.fasta(fasta_file, as.string = TRUE, forceDNAtolower=FALSE)
    cell4a <- readDNAStringSet(fasta_file, "fasta")
    names(cell4a) <- names(names.cell4)
   
    #List of lists of sequence of reads that are assigned to the same cluster_id
    id_seq <- vector("list", length(subgroups))
    for (i in 1:length(subgroups)){
        reads <- id_reads[[i]]
        id_seq[[i]] <- cell4a[reads]
    }
    
    return(id_seq)
}
