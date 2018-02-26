#' @export
umiGroup <- function(UMI1, threshold, UMI2 = NULL, groups=NULL, flip=NULL) 
# Computes the distances between UMIs of different reads 
# and assigns them to the same group if they are low-distance.
# 
# written by Florian Bieberich
# with modifications by Aaron Lun
# created 24 November 2017
{
    Nreads <- length(UMI1)
    if (missing(groups)) {
        groups <- rep(1L, Nreads)
    }
    if (length(groups)!=Nreads) { 
        stop("lengths of 'groups' and 'UMI1' should be the same")
    }
    if (!is.null(UMI2) && Nreads!=length(UMI2)) {
        stop("lengths of 'UMI1' and 'UMI2' should be the same")
    }
    if (!is.null(flip) && Nreads!=length(flip)) {
        stop("lengths of 'flip' and 'UMI1' should be the same")
    }

    by.group <- split(seq_len(Nreads), groups)
    output <- integer(Nreads)
    last <- 0L
    cur.UMI2 <- NULL
    
    filt.group <- which(lapply(by.group, length)>1)
    for (f in filt.group) { 
        current <- by.group[[f]]
        cur.UMI1 <- UMI1[current]

        if (!is.null(UMI2)) {
            cur.UMI2 <- UMI2[current]
            if (!is.null(flip)) { 
                cur.UMI2 <- UMI2[current]
                cur.flip <- flip[current]
                cur.tmp <- cur.UMI1[cur.flip]
                cur.UMI1[cur.flip] <- cur.UMI2[cur.flip]
                cur.UMI2[cur.flip] <- cur.tmp
            }  
        }

        UMI.id <- .umiClust(UMI1=cur.UMI1, UMI2=cur.UMI2, threshold=threshold)
        UMI.id <- UMI.id + last
        output[current] <- UMI.id 
        last <- max(UMI.id)
    }  
    return(output) 
}

#' @importFrom Biostrings stringDist
#' @importFrom igraph make_graph membership components
.umiClust <- function(UMI1, UMI2, threshold) {
    # stringDist defaults to method="quality" for QualityScaledXStringSet inputs.
    # Otherwise it remains as method="levenstein".
    lev.dist <- as.vector(stringDist(UMI1)) 
    if (!is.null(UMI2)) {
        lev.dist <- lev.dist + as.vector(stringDist(UMI2))
    }
   
    # Identifying links between reads, based on combined distance below the threshold.
    to.link <- which(lev.dist <= threshold)
    if (length(to.link)==0L) { 
        return(seq_along(UMI1))
    }

    # Getting indices without expanding the dist object. 
    N <- length(UMI1)
    col.ends <- cumsum((N-1L):1L)
    current.col <- findInterval(to.link, col.ends, left.open=TRUE) + 1L
    row.starts <- seq_len(N-1L)
    current.row <- row.starts[current.col] + to.link - c(0L, col.ends)[current.col]  

    # Identifying the graph components.
    g <- make_graph(rbind(current.col, current.row), n=N)
    return(membership(components(g)))
}

#' @export
#' @importFrom BiocGenerics order
umiGroup2 <- function(UMI1, max.lev1 = 3, UMI2 = NULL, max.lev2 = max.lev1, min.qual=10) 
# Uses a greedy approach to define UMI groupings, using only the Levenshtein distances.
# Some acknowledgement of the quality of the base calls is obtained with 'min.qual'.
# 
# written by Aaron Lun
# created 25 February 2018
{
    UMI1 <- .safe_masker(UMI1, threshold=min.qual)
    o1 <- order(UMI1)
    out1 <- .Call(cxx_umi_group, UMI1, o1 - 1L, max.lev1)

    # Repeating for the second UMI, if it is available.
    if (!is.null(UMI2)) { 
        UMI2 <- .safe_masker(UMI2, threshold=min.qual)
        o2 <- order(UMI2)
        out2 <- .Call(cxx_umi_group, UMI2, o2 - 1L, max.lev2)
    }

    out1 <- out    
    all.lengths <- lengths(out1)
    lost <- logical(length(all.lengths))
    group <- integer(length(all.lengths))
    counter <- 1L
    chosen <- which.max(all.lengths)

    repeat {
        # Deciding what to perform pairwise alignments to.
        targets <- out1[[chosen]]
        is.lost <- lost[targets]
        
        if (any(is.lost)) { 
            targets <- targets[!is.lost]
            out1[[chosen]] <- targets
            all.lengths[chosen] <- length(targets)
            new.chosen <- which.max(all.lengths)
            if (all.lengths[new.chosen] > length(targets)) {
                chosen <- new.chosen
                next 
            }
        }

        # Performing those pairwise alignments.
#        cur.scores <- pairwiseAlignment(umis[targets], umis[chosen], gapOpening = gapOpening, gapExtension = gapExtension, # ..., 
#                                        scoreOnly=TRUE)
#        success <- cur.scores >= min.score1

        group[targets] <- counter
        counter <- counter + 1L  

        lost[targets] <- TRUE
        all.lengths[targets] <- 0L
#        print(c(length(targets), sum(success), sum(all.lengths==0L)))
        print(sum(all.lengths==0))
        if (all(all.lengths==0L)) {
            break
        }

#        # Eradicating the current choice from the failures.
#        failures <- targets[!success]
#        for (failed in failures) {
#            out1[[failed]] <- setdiff(out1[[failed]], chosen)
#        }
#        all.lengths[failures] <- all.lengths[failures] - 1L

        chosen <- which.max(all.lengths)
    }
    
    return(group)
}
