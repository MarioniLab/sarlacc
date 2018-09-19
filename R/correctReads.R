#' @export
#' @importFrom BiocParallel bplapply SerialParam
#' @importClassesFrom IRanges IntegerList
#' @importFrom S4Vectors DataFrame
#' @importFrom methods as
correctReads <- function(reads, UMI1, UMI2=NULL, mm.cmd="minimap2", mm.args = NULL, working.dir=NULL, min.identity=0.7,
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
        working.dir <- tempfile(".")
        dir.create(working.dir)
        on.exit(unlink(working.dir, recursive = TRUE))
    }

    if(!is.null(UMI2)){
        names(UMI2.copy) <- seq_along(reads)
    }
  
    cluster.list <- groupReads(reads, mm.cmd=mm.cmd, mm.args=mm.args, working.dir=working.dir, min.identity=min.identity)
    umi.subgroups <- bplapply(cluster.list, FUN = .umi_group, UMI1 = UMI1, UMI2 = UMI2, umi.args = group.args, BPPARAM = BPPARAM)
    subclustered <- unlist(mapply(FUN = split, x = cluster.list, f = umi.subgroups, SIMPLIFY = FALSE), recursive = FALSE) 
        
    # Performing MSA and constructing consensus sequences.
    aligned.reads <- do.call(multiReadAlign, c(list(reads, groups = subclustered, BPPARAM = BPPARAM), mra.read.args))
    read.copy <- do.call(consensusReadSeq, c(list(aligned.reads, BPPARAM = BPPARAM), cons.read.args))
    aligned.UMI1 <- do.call(multiReadAlign, c(list(UMI1, groups = subclustered, BPPARAM = BPPARAM), mra.umi1.args))
    UMI1.copy <- do.call(consensusReadSeq, c(list(aligned.UMI1, BPPARAM = BPPARAM), cons.umi1.args))
 
    if (!is.null(UMI2)) {
        aligned.UMI2 <- do.call(multiReadAlign, c(list(UMI2, groups = subclustered, BPPARAM = BPPARAM), mra.umi2.args))
        UMI2.copy <- do.call(consensusReadSeq, c(list(aligned.UMI2, BPPARAM = BPPARAM), cons.umi2.args))
    }
 
    # Reassigning the read ID of the first read in each group as read ID of final consensus sequence
    first <- vapply(subclustered, FUN = "[", i = 1, FUN.VALUE = 0L)
    names(UMI1.copy) <- names(read.copy) <- names(origins) <- names(reads)[first] 
      
    output <- DataFrame(reads = read.copy, UMI1 = UMI1.copy)
    if (!is.null(UMI2)) {
        UMI2.copy <- names(UMI1.copy)
        output$UMI2 <-  UMI2.copy
    }
    
    output$origins <- as(origins, "IntegerList")
    return(output)
}

.umi_group <- function(idx, UMI1, UMI2, umi.args) {
    do.call(umiGroup, c(UMI1=UMI1[idx], UMI2=UMI2[idx], umi.args))
}
