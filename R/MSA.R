MSA <- function(id_seq_len.fil)

#Returns a DNAStringSet of multiple sequence alignments for the sequences from each cluster id

{
    msalign <- vector("list", length(id_seq_len.fil))
    
    for (i in 1:length(id_seq_len.fil)){
        msalign[[i]] <- muscle(id_seq_len.fil[[i]])
        msalign[[i]] <- DNAStringSet(msalign[[i]])
    }
    
    return(msalign)
}
