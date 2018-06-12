#' @export
#' @importFrom Biostrings writeQualityScaledXStringSet
#' @importFrom BiocParallel SerialParam
minimapMerge <- function(reads, UMI1, UMI2=NULL, mm.cmd="minimap2", mm.args = NULL, working.dir=NULL, min.identity=0.7, max.iter=3L,
        mra.read.args=list(), mra.umi1.args=mra.read.args, mra.umi2.args=mra.read.args, 
        cons.read.args=list(), cons.umi1.args=cons.read.args, cons.umi2.args=cons.read.args,
        group.args=list(), BPPARAM=SerialParam()) 
# Merges reads into final consensus sequences using pairwise alignments reported by minimap2.
#
# written by Cheuk-Ting Law
# with modifications by Aaron Lun
# created 12 June 2018
{
    if (!is.null(working.dir)) {
        working.dir <- tempfile()
        dir.create(working.dir)
        on.exit(unlink(working.dir, recursive=TRUE))
    }
   
    fpath <- file.path(working.dir, "reads.fastq")
    all.args <- c(mm.args, "-x ava-ont", "-c", fpath, fpath)
    nreads <- rep(1L, length(reads))

    iterations <- 0L
    while (iterations <= max.iter) {
        # Aligning with minimap2 to define read clusters.
        writeQualityScaledXStringSet(reads, fpath)
        raw.paf <- system2(mm.cmd, args = all.args, stdout = TRUE)
        cleaned.paf <- .process_paf(paf.raw, min.match=min.identity)
        cluster.list <- .cluster_paf(cleaned.paf, names(reads))
       
        # Defining UMI subclusters with umiGroup2.
        umi.subgroups <- bplapply(cluster.list, FUN=.umi_group, UMI1=UMI1, UMI2=UMI2, umi.args=group.args, BPPARAM=BPPARAM)
        subclustered <- unlist(mapply(FUN=split, x=cluster.list, f=umi.subgroups, SIMPLIFY=FALSE), recursive=FALSE)

        # Using the name of the first read as the name for each group.
        names(subclustered) <- names(reads)[vapply(cluster.list, FUN="[", i=1, FUN.VALUE=0L)]

        # Creating alignments and consensus sequences, and overwriting existing values.
        aligned.reads <- do.call(multiReadAlign, c(list(reads, groups=subclustered, BPPARAM=BPPARAM), mra.read.args))
        reads <- do.call(consensusReadSeq, c(list(aligned.reads, BPPARAM=BPPARAM), cons.read.args))
        
        aligned.UMI1 <- do.call(multiReadAlign, c(list(UMI1, groups=subclustered, BPPARAM=BPPARAM), mra.umi1.args))
        UMI1 <- do.call(consensusReadSeq, c(list(aligned.UMI1, BPPARAM=BPPARAM), cons.umi1.args))

        if (!is.null(UMI2)) {
            aligned.UMI2 <- do.call(multiReadAlign, c(list(UMI2, groups=subclustered, BPPARAM=BPPARAM), mra.umi2.args))
            UMI2 <- do.call(consensusReadSeq, c(list(aligned.UMI2, BPPARAM=BPPARAM), cons.umi2.args))
        }

        # Deciding whether to terminate or not.
        iterations <- iterations + 1L
        if (all(lengths(subclustered)==1L)) { 
            break
        }
        nreads <- unlist(lapply(subclustered, FUN=function(idx) { sum(nreads[idx]) }))
    }

    # Returning a list of values as output.
    output <- list(reads=aligned, UMI1=UMI1)
    if (!is.null(UMI2)) { 
        output$UMI2 <- UMI2
    } 
    output$nreads <- nreads
    return(output)
}

#' @importFrom igraph make_graph components V
.cluster_paf <- function(paf, all.names) {
    edges <- rbind(match(paf$qname, all.names), match(paf$tname, all.names))
    G <- make_graph(edges)

    comp <- components(G)$membership
    vertices <- as.vector(V(G))

    c(unname(split(vertices, comp)), 
      as.list(seq_along(all.names)[-vertices]))
}

#' @importFrom data.table fread
.process_paf <- function(paf, min.match=0.7) {
    paf <- fread(paf, select = c(1, 2, 6, 7, 10), sep="\t", header= FALSE)
    colnames(paf) <- c("qname", "qlength", "tname", "tlength", "alength")

    mean.length <- (paf$qlength+paf$tlength)/2
    prop.identical <- paf$matches/mean.length
    paf$identity <- prop.identical

    paf[prop.identical >= min.match,]
}

.umi_group <- function(idx, UMI1, UMI2, umi.args) {
    do.call(umiGroup2, c(UMI1=UMI1[idx], UMI2=UMI2[idx], umi.args))
}
