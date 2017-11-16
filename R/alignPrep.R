#The read data already needs to be aligned to the reference genome.
#therefore the bash script stored in /init has to be executed before running alignPrep.
#In order to be then read into the pipeline, .sam file has to be converted into .txt that only contains the first six columns:
#bahs command: "awk '{print $1,$2,$3,$4,$5,$6}' alignment.sam > align_pos.txt" #first 6 columns in new file.
alignPrep <- function(SAM_file)

#Returns an GRanges object containing read name, position, length, strand and chr location  

{
  
    alignment = read.table(SAM_file, skip=1)
    #table modification
    colnames(alignment) <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "SEQ")
    Read_len <- cigarWidthAlongReferenceSpace(alignment[,6]) #extract length from CIGAR string without softclips
    Query_len <- cigarWidthAlongQuerySpace(alignment[,6]) #with softclipped regions
    align_len <- cbind(alignment, Read_len, Query_len) #append to dataframe
    align_len$POS <- as.numeric(as.character(align_len$POS))
    align_len <- align_len[complete.cases(align_len[ ,4]),] #removes NA from the 4th column
    align_len[,3] <- gsub("^", "chr", align_len[,3])#extend RNAME with "chr" to fit GRanges format
    
    #make GRanges objects
    grange = with(align_len, GRanges(RNAME, IRanges(start = POS, width = Read_len, names = QNAME), strand = ifelse(FLAG=="0", "+", "-")))   
    
    
    return(list(grange = grange, alignment_data = align_len))
}


