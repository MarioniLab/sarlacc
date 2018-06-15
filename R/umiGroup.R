#' @export
#' @importFrom BiocGenerics order
umiGroup <- function(UMI1, max.lev1 = 3, UMI2 = NULL, max.lev2 = max.lev1)
# Groups UMIs based on their Levenshtein distances.
# Note that the 'importFrom' is due to an implicit order() call in the C++ code.
# 
# written by Aaron Lun
# created 25 February 2018
{
    out1 <- .Call(cxx_umi_group, UMI1, max.lev1)

    # Repeating for the second UMI, if it is available.
    if (!is.null(UMI2)) { 
        out2 <- .Call(cxx_umi_group, UMI2, max.lev2)
        out1 <- mapply(intersect, out1, out2)
    }

    .Call(cxx_descending_graph_cluster, out1) 
}



