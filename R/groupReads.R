#' @export
#' @importFrom Biostrings writeQualityScaledXStringSet subseq
groupReads <- function(reads, truncate=1000, working.dir=NULL, mm.cmd="minimap2", mm.args=NULL, min.identity=0.7)
# Uses minimap2 to perform pairwise alignment of reads.
{
    if (is.null(working.dir)) {
        working.dir <- tempfile(".")
        dir.create(working.dir, showWarnings=FALSE)
        on.exit({ unlink(working.dir, recursive=TRUE) })
    }

    fpath <- file.path(working.dir, "reads.fastq")
    all.args <- c(mm.args, "-x ava-ont", "-c", fpath, fpath)
    paf.cmd <- paste(c(mm.cmd, all.args, "|cut -f 1-12"), collapse = " ")

    if (is.finite(truncate)) {
        reads <- subseq(reads, 1L, width=pmin(width(reads), truncate))
    }
    writeQualityScaledXStringSet(reads, fpath)
    cleaned.paf <- .process_paf(paf.cmd, min.match = min.identity) # supply the minimap2 shell command directly to fread().
    .cluster_paf(cleaned.paf, sub(" .*", "", names(reads)))
}

#' @importFrom igraph make_graph components V
.cluster_paf <- function(paf, all.names) {
    edges <- rbind(match(paf$qname, all.names), match(paf$tname, all.names))
    if (length(edges)==0L) {
        return(as.list(seq_along(all.names)))
    }

    G <- make_graph(edges)
    comp <- components(G)$membership
    vertices <- as.vector(V(G))
    c(unname(split(vertices, comp)), as.list(seq_along(all.names)[-vertices]))
}

#' @importFrom data.table fread
.process_paf <- function(paf, min.match=0.7) {
    paf <- fread(paf, select = c(1, 2, 6, 7, 10), sep="\t", header= FALSE,fill=TRUE)
    colnames(paf) <- c("qname", "qlength", "tname", "tlength", "matches")
    mean.length <- (paf$qlength+paf$tlength)/2
    prop.identical <- paf$matches/mean.length
    paf$identity <- prop.identical
    paf[prop.identical >= min.match,]
}
