MSA <- function(cluster.reads)

#Returns a DNAStringSet of multiple sequence alignments for the sequences from each cluster id

{
    msalign <- vector("list", length(cluster.reads))
    
    for (i in 1:length(cluster.reads)){
        msalign[[i]] <- muscle(cluster.reads[[i]])
        msalign[[i]] <- DNAStringSet(msalign[[i]])
    }
    
    return(msalign)
}
