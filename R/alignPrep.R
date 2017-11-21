#The read data already needs to be aligned to the reference genome.
#therefore the bash script stored in /init has to be executed before running alignPrep.
#In order to be then read into the pipeline, .sam file has to be converted into .txt that only contains the first six columns:
#bash command: "awk '{print $1,$2,$3,$4,$5,$6,$10}' alignment.sam > alignment.txt" #first 6 columns in new file.

alignPrep <- function(sam, minq = 10)
#Returns an GRanges object containing read name, position, length, strand and chr location  
{
    # Figuring out how many headers to skip.
    N <- 0L
    curfile <- file(sam, open="r")
    repeat {
        curline <- readLines(curfile, n=1)
        if (!grepl("^@", curline)) {
            break
        }
        N <- N + 1L
    }
    close(curfile)

    # Setting "100" to ignore everything afterwards.
    what <- c(list("character","integer","character", "character", "integer", "character"), vector("list", 100))
    suppressWarnings(mapping <- read.table(sam, skip=N, colClasses=what, fill=TRUE, stringsAsFactors=FALSE))
    colnames(mapping) <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR")
   
    # Keeping only mapped reads and non-secondary reads. 
    keep <- !bitwAnd(mapping$FLAG, 0x4) & !bitwAnd(mapping$FLAG, 0x100)
    if (!is.null(minq)) {
        keep <- keep & mapping$MAPQ >= minq
    }
    mapping <- mapping[keep,]

    # Creating a GRanges object.
    align.len <- cigarWidthAlongReferenceSpace(mapping$CIGAR)
    pos <- as.integer(as.character(mapping$POS))
    granges <- GRanges(mapping$RNAME, IRanges(pos, width=align.len), strand=ifelse(bitwAnd(mapping$FLAG, 0x10), "-", "+"))

    # Splitting by the mapping names
    return(split(granges, mapping$QNAME))    
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
