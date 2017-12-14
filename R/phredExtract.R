# Extract phredscore of UMI

phredExtract <- function(align.stats, quality, adaptor2=NULL, position=NULL, length=NULL) 
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
    if(!length == position[2]-position[1]+1){
        stop("Please enter correct length")
    }
    length <- NULL
    
    # Identifying the number of deletions in the pattern before the start of the game.
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
    
    all.dels <- gregexpr("-", align.stats$pattern)
    indel.start <- indel.end <- integer(nrow(align.stats))
    for (i in seq_len(nrow(align.stats))) {
        dels <- as.integer(all.dels[[i]])
        if (length(dels)==1L && dels==-1L) { 
            next
        }
        true.pos <- dels - seq_along(dels)
        indel.start[i] <- sum(true.pos < position[1])
        indel.end[i] <- sum(true.pos < position[2])
    }
   
    if(!is.null(adaptor2)){
        cur.start <- align.stats$start.pattern
        align.stats$start.pattern <- align.stats$end.pattern
        align.stats$end.pattern <- align.stats$start.pattern
    }
    
    qual <- substr(aligned$quality, start=align.stats$start.pattern, stop=align.stats$end.pattern)
    
    if(!is.null(adaptor2)){
        qual <- reverse(qual)
    }
    
    phred.qual <- substr(qual, start=bump.start+position[1]-indel.start, stop=bump.end+position[2]-indel.end)
    
    return(PhredQuality(phred.qual))
    
}

# Get mean value for every phred - UMI string
#phred.mean <- numeric(length(extr))
#for (i in 1:length(extr)){
#    avg <- mean(strtoi(charToRaw(extr[i]),16L))-33
#    phred.mean[i] <- 10^(-avg/10)
#}


