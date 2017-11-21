#The read data already needs to be aligned to the reference genome.
#therefore the bash script stored in /init has to be executed before running alignPrep.
#In order to be then read into the pipeline, .sam file has to be converted into .txt that only contains the first six columns:
#bash command: "awk '{print $1,$2,$3,$4,$5,$6,$10}' alignment.sam > alignment.txt" #first 6 columns in new file.
alignPrep <- function(SAM_as_txt_file, txdb = NULL)

#Returns an GRanges object containing read name, position, length, strand and chr location  

{
  
    mapping = read.table(SAM_as_txt_file, skip=1)
    
    # Table modification.
    colnames(mapping) <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "SEQ")
    Read_len <- cigarWidthAlongReferenceSpace(mapping[,6]) # Extract length from CIGAR string without softclips.
    Query_len <- cigarWidthAlongQuerySpace(mapping[,6]) # With softclipped regions.
    align_len <- cbind(mapping, Read_len, Query_len)
    align_len$POS <- as.numeric(as.character(align_len$POS))
    align_len <- align_len[complete.cases(align_len[ ,4]),] # Removes NA from the 4th column.
    align_len[,3] <- gsub("^", "chr", align_len[,3])# Extend RNAME with "chr" to fit GRanges format.
    
    # Transform table to GRanges object.
    grange = with(align_len, GRanges(RNAME, IRanges(start = POS, width = Read_len, names = QNAME), strand = ifelse(FLAG=="0", "+", "-")))   
    
    if(!is.null(txdb)){
        grange <- .mColAdd(grange = grange, txdb = txdb)
    }
    
    return(list(grange = grange, mapping_data = align_len))
}


.mColAdd <- function(grange, txdb)
    # Adds gene_id containing metacolumn to grange.
{
    exon_data <- exons(txdb)
    gene_id <- mapIds(txdb, keytype="EXONID", key= as.character(exon_data$exon_id), column="GENEID")
    
    # Overlap between mapped reads and annotated exons.
    alignment_gr_reads <- findOverlaps(grange, exon_data)
    
    # Assign gene ids as metacolumn to gr4 and add NA for reads without an identified gene id
    
    gene_id_uni <- gene_id[subjectHits(alignment_gr_reads[which(duplicated(queryHits(alignment_gr_reads))==FALSE)])]
    
    gene_id_allread <- rep(NA, length(grange))
    gene_id_allread[unique(queryHits(alignment_gr_reads))] <- gene_id_uni
    mcols(grange) <- gene_id_allread
    names(mcols(grange)) <- "gene_id"
    gr_gene <- grange
    
    return(gr_gene)
}
