#The read data already needs to be aligned to the reference genome.
#therefore the bash script stored in /init has to be executed before running alignPrep.

alignPrep <- function(SAM_file)

#Returns an GRanges object containing read name, position, length, strand and chr location  

{
  
    alignments_cell4 = read.table(SAM_file, skip=1)
    #table modification
    colnames(alignments_cell4) <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR")
    Read_len <- cigarWidthAlongReferenceSpace(alignments_cell4[,6]) #extract length from CIGAR string without softclips
    Query_len <- cigarWidthAlongQuerySpace(alignments_cell4[,6]) #with softclipped regions
    align_cell4_len <- cbind(alignments_cell4, Read_len, Query_len) #append to dataframe
    align_cell4_len$POS <- as.numeric(as.character(align_cell4_len$POS))
    align_cell4_len <- align_cell4_len[complete.cases(align_cell4_len[ ,4]),] #removes NA from the 4th column
    align_cell4_len[,3] <- gsub("^", "chr", align_cell4_len[,3])#extend RNAME with "chr" to fit GRanges format
    
    #make GRanges objects
    grange = with(align_cell4_len, GRanges(RNAME, IRanges(start = POS, width = Read_len, names = QNAME), strand = ifelse(FLAG=="0", "+", "-")))   
    
    
    return(grange)
}
