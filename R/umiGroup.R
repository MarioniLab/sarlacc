#' @export
#' @importFrom BiocGenerics order
#' @importFrom igraph make_graph components
umiGroup <- function(UMI1, max.lev1 = 3, UMI2 = NULL, max.lev2 = max.lev1, max.err=NA, use.densities=TRUE)
# Groups UMIs based on their Levenshtein distances.
# Note that the 'importFrom' is due to an implicit order() call in the C++ code.
# 
# written by Aaron Lun
# created 25 February 2018
{
    if (length(UMI1)==1L) { 
        return(1)
    }

    UMI1 <- .safe_masker(UMI1, max.err)
    out1 <- .Call(cxx_umi_group, UMI1, max.lev1)

    # Repeating for the second UMI, if it is available.
    if (!is.null(UMI2)) {
        if (length(UMI1)!=length(UMI2)) { 
            stop("mismatch in lengths between 'UMI1' and 'UMI2'")
        }
        UMI2 <- .safe_masker(UMI2, max.err)
        out2 <- .Call(cxx_umi_group, UMI2, max.lev2)
        out1 <- mapply(intersect, out1, out2, SIMPLIFY=FALSE)
    }

    # if (use.densities) {
    #     return(.Call(cxx_descending_graph_cluster, out1))
    # } else {
    #     left <- unlist(out1) + 1L
    #     right <- rep(seq_along(out1), lengths(out1))
    #     G <- make_graph(rbind(left, right), directed=FALSE, n=length(UMI1))
    #     return(components(G)$membership)
    # }
    .Call(cxx_cluster_umis, out1)
}

#' @importFrom Biostrings quality
#' @importFrom methods is as
#' @importClassesFrom Biostrings DNAStringSet QualityScaledDNAStringSet
#' @importClassesFrom IRanges NumericList
.safe_masker <- function(UMI, max.err) {
    has.quals <- is(UMI, "QualityScaledDNAStringSet")
    if (!is.na(max.err) && has.quals) { 
        all.qual <- as.list(as(quality(UMI), "NumericList"))
        UMI <- .Call(cxx_mask_bad_bases, UMI, all.qual, max.err)
    } else if (has.quals) {
        UMI <- as(UMI, "DNAStringSet")
    } 
    return(UMI)
}

.central_max <- function(out){
    out.groups <- rep(seq_along(out), lengths(out))
    flat.out <- unlist(out)
    collected.group <- c()
    collected.member <- c()
    
    for(i in seq_along(out)){
        member.count <- rle(sort(out.groups))
        member.num <- member.count$lengths
        member <- member.count$values
        max.group <- member[which.max(member.num)]
        target <- out.groups==max.group
        if(length(target)==0){
            break
        }
        collected.member <- c(collected.member, flat.out[target])
        collected.group <- c(collected.group, out.groups[target])
        remove <- flat.out%in%flat.out[target]
        flat.out <- flat.out[!remove]
        out.groups <- out.groups[!remove]
    }
    group.order <- order(collected.member)
    collected.group[group.order]
}
