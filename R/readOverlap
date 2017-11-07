readOverlap <- function(grange, minoverlap = 100)

#Returns overlapped reads as query-subject overlap and as splitted list for query reads

{
    #Comparison of overlap
    overlap <- findOverlaps(grange, minoverlap = minoverlap)
    splitted_overlap <- split(subjectHits(overlap),queryHits(overlap))
    
    return(list(overlap, splitted_overlap))
}
