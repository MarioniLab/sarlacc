#' @export
umiGroup <- function(UMI1, threshold1 = 3, UMI2 = NULL, threshold2 = threshold1, max.err=NA, groups=NULL)
# Groups UMIs based on their Levenshtein distances.
# Note that the 'importFrom' is due to an implicit order() call in the C++ code.
# 
# written by Aaron Lun
# created 25 February 2018
{
    UMI1 <- qualityMask(UMI1, max.err)
    if (!is.null(UMI2)) {
        UMI2 <- qualityMask(UMI2, max.err)
    }

    if (is.null(groups)) { 
        by.group <- list(seq_along(UMI1))
    } else if (!is.list(groups)) {
        by.group <- split(seq_along(UMI1), groups)
    } else {
        by.group <- groups
    }

    out <- .Call(cxx_umi_group, UMI1, threshold1, UMI2, threshold2, by.group)
    unlist(out, recursive=FALSE, use.names=FALSE)
}
