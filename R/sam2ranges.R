#' @export
#' @importFrom utils count.fields read.delim 
#' @importFrom GenomicAlignments cigarWidthAlongReferenceSpace
#' @importFrom S4Vectors split
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb Seqinfo
sam2ranges <- function(sam, minq = 10, restricted = NULL)
# Returns an GRanges object containing read name, position, length, strand and chr location.
# This is a bit more complicated due to the need to parse a SAM file, not a BAM file 
# (which would have been easy via Rsamtools, but ONT CIGARs don't fit inside BAM fields).
#
# written by Florian Bieberich
# with modifications by Aaron Lun
{
    # Figuring out how many headers to skip.
    N <- 1L
    curfile <- file(sam, open="r")
    collected <- c()
    is.empty <- FALSE
    repeat {
        curline <- readLines(curfile, n=1)
        if (length(curline)==0L) {
            is.empty <- TRUE
            break
        }
        if (!grepl("^@", curline)) {
            break
        }
        if (grepl("^@SQ", curline)) {
            collected[N] <- curline
        }
        N <- N + 1L
    }
    close(curfile)

    ref.len <- as.integer(sub(".*\tLN:([^\t]+)(\t.*)?","\\1", collected))
    names(ref.len) <- sub(".*\tSN:([^\t]+)(\t.*)?", "\\1", collected)
    ref.len <- c(ref.len, `*`=0L) # for unmapped reads.

    # Avoid silly behaviour if the file is actually empty.
    if (is.empty) { 
        out <- GRangesList()
        seqinfo(out) <- Seqinfo(names(ref.len), seqlengths=ref.len)
        return(out)
    }

    # Actually extracting out the data.
    nfields <- max(count.fields(sam, skip=N, sep="\t", comment.char="", quote=""), na.rm=TRUE)
    what <- list("character","integer","character", "character", "integer", "character")
    what <- c(what, vector("list", nfields - length(what)))
    suppressWarnings(mapping <- read.delim(sam, header=FALSE, skip=N, colClasses=what, fill=TRUE, quote="", comment.char="", stringsAsFactors=FALSE))
    colnames(mapping) <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR")
   
    # Only processing mapped reads above a certain quality and on the desired chromosomes.
    keep <- !bitwAnd(mapping$FLAG, 0x4) 
    if (!is.null(minq)) {
        keep <- keep & mapping$MAPQ >= minq
    }
    if (!is.null(restricted)){ 
        keep <- keep & mapping$RNAME %in% restricted
    }

    mapping$RNAME[!keep] <- "*"
    mapping$POS[!keep] <- "1" # avoid trim() warnings.
    mapping$CIGAR[!keep] <- ""

    strandedness <- ifelse(bitwAnd(mapping$FLAG, 0x10), "-", "+")
    strandedness[!keep] <- "*"

    # Creating a GRanges object.
    align.len <- cigarWidthAlongReferenceSpace(mapping$CIGAR)
    pos <- as.integer(mapping$POS)
    granges <- GRanges(mapping$RNAME, IRanges(pos, width=align.len), strand=strandedness,
        seqinfo=Seqinfo(names(ref.len), seqlengths=ref.len))

    granges$left.clip <- .get_clip_length(mapping$CIGAR)
    granges$right.clip <- .get_clip_length(mapping$CIGAR, start=FALSE)

    split(granges, mapping$QNAME, drop=FALSE)
}

.get_clip_length <- function(cigars, start=TRUE) {
    cliplen <- integer(length(cigars))
    for (op in c("H", "S")) { # hard clips before soft clips.
        if (start) {
            finder <- paste0("^[0-9]+", op)
            keeper <- paste0("^([0-9]+)", op, ".*")
        } else {
            finder <- paste0("[0-9]+", op, "$")
            keeper <- paste0(".*[^0-9]([0-9]+)", op, "$")
        }
        present <- grepl(finder, cigars)
        cliplen[present] <- cliplen[present] + as.integer(sub(keeper, "\\1", cigars[present]))
        cigars <- sub(finder, "", cigars)
    }
    return(cliplen)
}

