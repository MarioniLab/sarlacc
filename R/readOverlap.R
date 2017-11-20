readOverlap <- function(map_data, minoverlap = 100)

#Returns overlapped reads as query-subject overlap and as splitted list for query reads

{
    #Comparison of overlap
    overlap <- findOverlaps(map_data$grange, minoverlap = minoverlap, ignore.strand = TRUE)
    splitted_overlap <- split(subjectHits(overlap),queryHits(overlap))
    
    return(list(overlap, splitted_overlap))
}
