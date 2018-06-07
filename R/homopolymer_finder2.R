homopolymer_finder <- function(seq_str){
    seq_letters <- unique(strsplit(seq_str, "", perl = TRUE)[[1]])
    seq_str <- DNAString(as.character(seq_str))
    homo_range <- bplapply(seq_letters,function(x){
        match_loc <- matchPattern(x, seq_str)
        repeat_group <- t(.consencutive_finder(match_loc@ranges@start))
        cbind.data.frame(repeat_group, x)
    })
    homo_range <- data.table(do.call(rbind.data.frame, homo_range))
    colnames(homo_range) <- c("start", "end", "letters")
    homo_ir <- IRanges(start = homo_range$start, end = homo_range$end)
    elementMetadata(homo_ir)$letters <- homo_range$letters
    homo_ir <- homo_ir[order(homo_ir@start)]
    homo_ir[homo_ir@width != 1 & homo_ir@elementMetadata$letters != "N"]
}


.consencutive_finder <- function(vec){
    breaks <- c(0, which(diff(vec) != 1), length(vec)) 
    sapply(seq(length(breaks) - 1), 
           function(x)vec[c((breaks[x] + 1), breaks[x+1])])
}

homopolymer_finder("ATGGGGATTGCAAAACCCCCAAAAAAAAAAA")
