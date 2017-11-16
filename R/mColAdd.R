mColAdd <- function(grange, txdb)

#Returns Granges object as initially made but with an additional meta column containing the corresponding gene id

{
    exon_data <- exons(txdb)
    gene_id <- mapIds(txdb, keytype="EXONID", key= as.character(exon_data$exon_id), column="GENEID")
  
    alignment_gr_reads <- findOverlaps(grange, exon_data) #exon id to respective read (alignment)
  
    #Assign gene ids as metacolumn to gr4 and add NA for reads without an identified gene id
    gene_id_uni <- gene_id[subjectHits(alignment_gr_reads[which(duplicated(queryHits(alignment_gr_reads))==FALSE)])]
  
    gene_id_allread <- vector("character", length(grange))
    gene_id_allread[gene_id_allread==""] <- NA
    gene_id_allread[unique(queryHits(alignment_gr_reads))] <- gene_id_uni
    mcols(grange) <- gene_id_allread
    names(mcols(grange)) <- "gene_id"
    gr_gene <- grange
    
    return(gr_gene)
}

