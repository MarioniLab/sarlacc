readOverlap <- function(map_data, minoverlap = 100)

#Returns overlapped reads as query-subject overlap and as splitted list for query reads

{
    #Comparison of overlap
    overlap <- findOverlaps(map_data, minoverlap = minoverlap, ignore.strand = TRUE)
    overlap.reads <- split(subjectHits(overlap),queryHits(overlap))
    
    return(overlap.reads)
}
