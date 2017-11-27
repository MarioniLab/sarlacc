# For use after findOverlaps on the read positions. 
# Filters reads for each subgroup.

clusterReads <- function(overlaps) {
    s <- subjectHits(overlaps)
    q <- queryHits(overlaps)
    keep <- s!=q
    g <- make_graph(rbind(s[keep], q[keep]), n=queryLength(overlaps))
    return(components(g)$membership)
}

