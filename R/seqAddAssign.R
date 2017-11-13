seqAddAssign <- function(fastq_file, subgroups, id_reads)

#Returns list of cluster id groups that contain the corresponding sequences

{
    fastq <- readFastq(fastq_file, withIds=TRUE)
    seqs <- sread(fastq)
    names(seqs) <- ShortRead::id(fastq)
   
    #List of lists of sequence of reads that are assigned to the same cluster_id
    id_seq <- vector("list", length(subgroups))
    for (i in 1:length(subgroups)){
        reads <- id_reads[[i]]
        id_seq[[i]] <- seqs[reads]
    }
    
    return(id_seq)
}

