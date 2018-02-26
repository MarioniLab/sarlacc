#' @export
#' @importFrom utils count.fields read.table
#' @importFrom GenomicAlignments cigarWidthAlongReferenceSpace
#' @importFrom S4Vectors split
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
alignPrep <- function(sam, minq = 10, restricted = NULL)
# Returns an GRanges object containing read name, position, length, strand and chr location.
# This is a bit more complicated due to the need to parse a SAM file, not a BAM file 
# (which would have been easy via Rsamtools, but ONT CIGARs don't fit inside BAM fields).
#
# written by Florian Bieberich
# with modifications by Aaron Lun
{
    # Figuring out how many headers to skip.
    N <- 0L
    curfile <- file(sam, open="r")
    repeat {
        curline <- readLines(curfile, n=1)
        if (!grepl("^@", curline)) {
            break
        }
        N <- N + 1L
    }
    close(curfile)

    # Setting "100" to ignore everything afterwards.
    nfields <- max(count.fields(sam, skip=N, sep="\t", comment.char="", quote=""), na.rm=TRUE)
    what <- list("character","integer","character", "character", "integer", "character")
    what <- c(what, vector("list", nfields - length(what)))
    suppressWarnings(mapping <- read.table(sam, skip=N, colClasses=what, fill=TRUE, quote="", comment.char="", stringsAsFactors=FALSE))
    colnames(mapping) <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR")
   
    # Keeping only mapped reads and non-secondary reads. 
    keep <- !bitwAnd(mapping$FLAG, 0x4) & !bitwAnd(mapping$FLAG, 0x100)
    if (!is.null(minq)) {
        keep <- keep & mapping$MAPQ >= minq
    }
    mapping <- mapping[keep,]

    # Creating a GRanges object.
    align.len <- cigarWidthAlongReferenceSpace(mapping$CIGAR)
    pos <- as.integer(as.character(mapping$POS))
    granges <- GRanges(mapping$RNAME, IRanges(pos, width=align.len), strand=ifelse(bitwAnd(mapping$FLAG, 0x10), "-", "+"))

    # Restricting to a subset of the alignments on particular chromosomes.
    read.names <- mapping$QNAME
    if (!is.null(restricted)){ 
        keep <- seqnames(granges) %in% restricted
        granges <- granges[keep]
        read.names <- factor(read.names)
        read.names <- read.names[as.logical(keep)]
    }
        
    split(granges, read.names, drop=FALSE)
}

