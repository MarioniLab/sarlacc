homopolymerMatcher <- function(aln, ref){
    homo_ir <- homopolymerFinder(as.character(ref))[[1]]
    aln_range <- IRanges(start = start(pattern(aln)), end = end(pattern(aln)))
    overlap_range <- findOverlaps(aln_range, homo_ir)
    homo_range <- lapply(split(overlap_range@to, f = overlap_range@from),range)
    homo_regions <- lapply(homo_range, function(x)1:length(homo_ir)%in%(x[1]:x[2]))
    homo_regions <- do.call(rbind,homo_regions)
    aligned_aln <- aligned(aln)

    homo_mat <- matrix(as.character(unlist(extractAt(aligned_aln, homo_ir))),,length(homo_ir),byrow = TRUE)
    homo_mat[!homo_regions] <- NA
    homo_mat <- apply(homo_mat, 2, function(x)do.call(rbind,strsplit(x, split = "", perl = TRUE)))
    homo_length <- sapply(homo_mat, ncol)
    match_list <- lapply(1:length(homo_length), function(x)homo_mat[[x]]==homo_ir@elementMetadata$letters[x])
    match_count <- lapply(match_list,rowSums)
    match_count <- lapply(1:length(match_count), function(x)factor(match_count[[x]], levels = 0:homo_length[x]))
    elementMetadata(homo_ir)$count <- NumericList(lapply(match_count, table))
    homo_ir
}
