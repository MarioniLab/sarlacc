#' @export
#' @importFrom IRanges findOverlaps coverage
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom BiocGenerics width strand<-
#' @importFrom methods is as
#' @importClassesFrom GenomicRanges GenomicRanges
#' @importFrom GenomeInfoDb seqnames
groupAlignedGenome <- function(locations) 
# Groups alignments based on their genomic coordinates.
#
# written by Aaron Lun
# created 22 November 2018
{
    mapped <- seqnames(locations)!="*"
    pcov <- coverage(locations[mapped & strand(locations)=="+"])
    ncov <- coverage(locations[mapped & strand(locations)=="-"])

    # Creating 'standard' intervals to overlap against.
    pstandards <- as(slice(pcov, lower=1), "GenomicRanges")
    nstandards <- as(slice(ncov, lower=1), "GenomicRanges")
    strand(pstandards) <- "*"
    strand(nstandards) <- "-"
    standards <- c(pstandards, nstandards)

    # All mapped elements should overlap with at least one standard,
    # and exactly one if the elements do not have multiple entries.
    olap <- findOverlaps(locations, standards, ignore.strand=FALSE, select="first")
    .clear_unmapped(olap, locations, mapped)
}

#' @export
#' @importFrom S4Vectors head
#' @importFrom BiocGenerics strand
#' @importFrom GenomeInfoDb seqnames
groupAlignedTrans <- function(locations) 
# Groups alignments based on the reference sequence name.
#
# written by Aaron Lun
# created 22 November 2018
{
    mapped <- seqnames(locations)!="*"
    first <- c(1L, 1L+cumsum(head(lengths(locations), -1)))
    by.rname <- unlist(seqnames(locations), use.names=FALSE)[first]
    by.strand <- unlist(strand(locations), use.names=FALSE)[first]
    .clear_unmapped(
        as.integer(factor(paste0(by.rname, ":", by.strand))), 
        locations, mapped)
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

