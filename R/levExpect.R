levExpect <- function(align.stats, x=100)
{
    align.stats <- align.stats[1:x,]
    all.dels <- gregexpr("-", align.stats$subject)
    bump.start <- bump.end <- integer(nrow(align.stats))
    for (i in seq_len(nrow(align.stats))) {
        dels <- as.integer(all.dels[[i]])
        if (length(dels)==1L && dels==-1L) { 
            next
        }
        true.pos <- dels - seq_along(dels)
        bump.start[i] <- sum(true.pos < position[1])
        bump.end[i] <- sum(true.pos < position[2])
    }

    a <- substr(align.stats$pattern, start=bump.start+position[1]-6, stop=bump.start+position[1]-1)
    b <- substr(align.stats$pattern, start=bump.end+position[2]+1, stop=bump.end+position[2]+6)
    
    a <- gsub("-", "", a)
    b <- gsub("-", "", b)

    curseq <- gsub("-", "", align.stats$subject[1])
    all.Ns <- gregexpr("N+", curseq)[[1]]
    position <- c(all.Ns, all.Ns+attr(all.Ns, "match.length")-1L)

    as <- substr(curseq, start=position[1]-6, stop=position[1]-1)
    bs <- substr(curseq, start=position[2]+1, stop=position[2]+6)

    lev.dist.a <- as.vector(stringDist(c(as,a), method = "levenshtein"))
    lev.dist.b <- as.vector(stringDist(c(bs,b), method = "levenshtein"))

    a.dist <- lev.dist.a[1:length(a)]
    b.dist <- lev.dist.b[1:length(b)]

    return(c(mean(a.dist+b.dist),sd(a.dist+b.dist)))
}
