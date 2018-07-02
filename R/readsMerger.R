#' @export
#' @importFrom Biostrings writeQualityScaledXStringSet
#' @importFrom BiocParallel SerialParam
#' @importClassesFrom IRanges IntegerList
#' @importFrom S4Vectors DataFrame
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
    if (is.null(working.dir)) {
        working.dir <- tempfile()
        dir.create(working.dir)
        on.exit(unlink(working.dir, recursive = TRUE))
    }
    fpath <- file.path(working.dir, "reads.fastq")
    all.args <- c(mm.args, "-x ava-ont", "-c", fpath, fpath)
    paf.cmd <- paste(c(mm.cmd, all.args, "|cut -f 1-12"), collapse = " ")
    origins <- as.list(seq_along(reads))
    read.copy <- reads
    UMI1.copy <- UMI1
    UMI2.copy <- UMI2
    names(origins) <- names(read.copy) <- names(UMI1.copy) <- seq_along(reads)
    
    if(!is.null(UMI2)){
    names(UMI2.copy) <- seq_along(reads)
    }
  
    iterations <- 0L
    while (iterations <= max.iter) {
        writeQualityScaledXStringSet(read.copy, fpath)
        cleaned.paf <- .process_paf(paf.cmd, min.match = min.identity)
        cluster.list <- .cluster_paf(cleaned.paf, names(read.copy))
        umi.subgroups <- bplapply(cluster.list, FUN = .umi_group, 
                                  UMI1 = UMI1.copy, UMI2 = UMI2.copy, umi.args = group.args, 
                                  BPPARAM = BPPARAM)
        subclustered <- unlist(mapply(FUN = split, x = cluster.list, 
                                      f = umi.subgroups, SIMPLIFY = FALSE), recursive = FALSE)
        
        previous.origins <- origins #Storing information from the previous round.
        previous.read.copy <- read.copy
        previous.UMI1.copy <- UMI1.copy
        origins <- lapply(subclustered, FUN = function(idx) {
          unlist(origins[idx], use.names = FALSE)
        })
    
        #Identifying the groups with changes and only performing MSA on those groups.
        names(origins) <- sapply(origins, FUN=function(x)paste(x, collapse = "_")) #Combining all origins as group ID
        names(previous.origins)<- sapply(previous.origins, FUN=function(x)paste(x, collapse = "_"))
        changed <- !names(origins)%in%names(previous.origins)
        changed.groups <- origins[changed]
        
        aligned.reads <- do.call(multiReadAlign, c(list(reads, groups = changed.groups, BPPARAM = BPPARAM), mra.read.args))
        read.copy <- do.call(consensusReadSeq, c(list(aligned.reads,BPPARAM = BPPARAM), cons.read.args))
        aligned.UMI1 <- do.call(multiReadAlign, c(list(UMI1, groups = changed.groups, BPPARAM = BPPARAM), mra.umi1.args))
        UMI1.copy <- do.call(consensusReadSeq, c(list(aligned.UMI1, BPPARAM = BPPARAM), cons.umi1.args))
        names(UMI1.copy) <- names(read.copy) <- names(changed.groups)
    
        #Combining the changed and unchanged groups.
        read.copy <- c(previous.read.copy[names(origins)[!changed]], read.copy)
        read.copy <- read.copy[order(names(read.copy))]
        UMI1.copy <- c(previous.UMI1.copy[names(origins)[!changed]], UMI1.copy)
        UMI1.copy <- UMI1.copy[order(names(UMI1.copy))]
        origins <- c(origins[!changed], origins[changed])
        origins <- origins[order(names(origins))]
        
        if (!is.null(UMI2)) {
          previous.UMI2.copy <- UMI2.copy
          aligned.UMI2 <- do.call(multiReadAlign, c(list(UMI2, groups = changed.groups, BPPARAM = BPPARAM), mra.umi2.args))
          UMI2.copy <- do.call(consensusReadSeq, c(list(aligned.UMI2, BPPARAM = BPPARAM), cons.umi2.args))
          names(UMI2.copy) <- names(changed.groups)
          UMI2.copy <- c(previous.UMI2.copy[names(origins)[!changed]], UMI2.copy)
          UMI2.copy <- UMI2.copy[order(names(UMI2.copy))]
        }
        
        iterations <- iterations + 1L
        if (all(lengths(subclustered) == 1L)) {
          break
        }
    }
  
      #Reassigning the read ID of the first read in each group as read ID
      names(UMI1.copy) <- names(read.copy) <- names(origins) <- names(reads)[vapply(origins, FUN = "[", i = 1, FUN.VALUE = 0L)] 
      
      output <- DataFrame(reads = read.copy, UMI1 = UMI1.copy)
      if (!is.null(UMI2)) {
        UMI2.copy <- names(UMI1.copy)
        output$UMI2 <-  UMI2.copy
      }
    
      output$origins <- as(origins, "IntegerList")
      return(output)
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

.umi_group <- function(idx, UMI1, UMI2, umi.args) {
    do.call(umiGroup, c(UMI1=UMI1[idx], UMI2=UMI2[idx], umi.args))
}
