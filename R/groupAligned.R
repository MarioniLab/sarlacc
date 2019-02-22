#' @export
#' @importFrom IRanges findOverlaps coverage slice
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom BiocGenerics width strand<- sort
#' @importFrom methods is as
#' @importClassesFrom GenomicRanges GenomicRangesList
#' @importFrom GenomeInfoDb seqnames
groupAlignedGenome <- function(locations) 
# Groups alignments based on their genomic coordinates.
#
# written by Aaron Lun
# created 22 November 2018
{
    if (is(locations, "GenomicRangesList")) {
        ref <- unlist(locations)
    } else {
        ref <- locations
    }
    
    mapped <- seqnames(ref)!="*"
    pstandards <- reduce(ref[mapped & strand(ref)=="+"])
    nstandards <- reduce(ref[mapped & strand(ref)=="-"])
    standards <- c(pstandards, nstandards)

    # All mapped elements should overlap with at least one standard,
    # and exactly one if the elements do not have multiple entries.
    if (is(locations, "GenomicRangesList")) {
        locations <- sort(locations)
    }
    olap <- findOverlaps(locations, standards, ignore.strand=FALSE, select="first")
    .clear_unmapped(olap, locations, mapped)
}

#' @export
#' @importFrom S4Vectors head
#' @importFrom BiocGenerics strand sort
#' @importFrom GenomeInfoDb seqnames
#' @importFrom methods is
#' @importClassesFrom GenomicRanges GenomicRangesList
groupAlignedTrans <- function(locations) 
# Groups alignments based on the reference sequence name.
#
# written by Aaron Lun
# created 22 November 2018
{
    if (is(locations, "GenomicRangesList")) {
        locations <- locations[seqnames(locations)!="*"]
        locations <- sort(locations)

        # Take the first entry within each GRL element 
        # (protect against out-of-range errors with lengths of zero).
        first <- cumsum(lengths(locations)) - lengths(locations) + (lengths(locations)>0L)
        by.rname <- unlist(seqnames(locations), use.names=FALSE)[first]
        by.strand <- unlist(strand(locations), use.names=FALSE)[first]
    } else {
        by.rname <- seqnames(locations)
        by.strand <- strand(locations)
    }
    
    combo <- paste0(as.integer(by.rname), ":", as.integer(by.strand))
    .clear_unmapped(as.integer(factor(combo)), locations, seqnames(locations)!="*")
}

#' @importFrom methods is
#' @importClassesFrom GenomicRanges GenomicRanges
.clear_unmapped <- function(grouping, locations, mapped) {
    if (is(locations, "GenomicRanges")) {
        grouping[-mapped] <- NA_integer_
    } else {
        grouping[!any(mapped)] <- NA_integer_
    }
    names(grouping) <- names(locations)
    grouping
}
