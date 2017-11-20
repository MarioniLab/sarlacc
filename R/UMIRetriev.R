UMIRetriev <- function(align.stats, position=NULL) 
# Pull out UMIs from the adaptor alignment data.
{
    if (is.null(position)) { 
        curseq <- gsub("-", "", align.stats$subject[1])
        all.Ns <- gregexpr("N+", curseq)[[1]]
        position <- c(all.Ns, all.Ns+attr(all.Ns, "match.length")-1L)
    } else {
        if (length(position)!=2L || position[1] > position[2] || position[1] <= 0L) {
            stop("invalid 'position' vector")
        }
    }
    
    # Identifying the number of deletions in the pattern before the start of the game.
    all.dels <- gregexpr("-", align.stats$subject)
    bump.start <- bump.end <- integer(nrow(align.stats))
    for (i in seq_len(nrow(align.stats))) {
        dels <- as.integer(all.dels[[i]])
        true.pos <- dels - seq_along(dels)
        bump.start[i] <- sum(true.pos < position[1])
        bump.end[i] <- sum(true.pos < position[2])
    }

    # Pulling out the UMI and cleaning it up.
    umi <- substr(align.stats$pattern, start=bump.start+position[1], stop=bump.end+position[2])
    umi <- gsub("-", "", umi)
    return(umi)
}


