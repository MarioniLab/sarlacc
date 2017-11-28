umiGroup <- function(UMI1, UMI2 = NULL, groups=NULL, threshold=2, flip=FALSE) 
# Function for both UMIs.
# umiGroup adds the levenshtein distances of both UMI's up and then proceeds as in the initial UMIgroup function.
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

    by.group <- split(seq_len(Nreads), groups)
    output <- integer(Nreads)
    last <- 0L

    for (f in names(by.group)) { 
        current <- by.group[[f]]

        if (is.null(UMI2) || !flip) { 
            UMI.id <- .umiClust(UMI1=UMI1[current], UMI2=UMI2[current], threshold=threshold)
        } else {
            UMI.id <- .umiCrossClust(UMI1=UMI1[current], UMI2=UMI2[current], threshold=threshold)
        }

        output[current] <- UMI.id + last
        last <- last + max(UMI.id)
    }  
    return(output) 
}

.umiClust <- function(UMI1, UMI2=NULL, threshold=2) {
    lev.dist <- as.vector(stringDist(UMI1, method = "levenshtein"))
    if (!is.null(UMI2)) {
        lev.dist <- lev.dist + as.vector(stringDist(UMI2, method = "levenshtein"))
    }
   
    # Identifying links between reads, based on combined distance below the threshold.
    to.link <- which(lev.dist <= threshold)
    if (length(to.link)==0L) { 
        return(seq_along(UMI1))
    }

    # Getting indices without expanding the dist structure.
    N <- length(UMI1)
    ind <- .getDistIndices(to.link, N)

    # Identifying the components of the graph.
    g <- make_graph(rbind(ind$col, ind$row), n=N)
    return(membership(components(g)))
}

.getDistIndices <- function(chosen, N) {
    col.ends <- cumsum((N-1L):1L)
    current.col <- findInterval(chosen, col.ends, left.open=TRUE) + 1L
    row.starts <- seq_len(N-1L)
    current.row <- row.starts[current.col] + chosen - c(0L, col.ends)[current.col]  
    return(list(row=current.row, col=current.col))
}

.umiCrossClust <- function(UMI1, UMI2, threshold=2) { 
    combined <- c(UMI1, UMI2)
    lev.dist <- as.vector(stringDist(combined, method="levenshtein"))
    N <- length(UMI1)

    # Figuring out the original entries.
    all.ind <- .getDistIndices(seq_along(lev.dist), length(combined))
    true.ind <- rep(seq_len(N), 2L)
    true.row <- true.ind[all.ind$row]
    true.col <- true.ind[all.ind$col]

    # Pulling out the parts of the matrix corresponding to a straight comparison
    # between UMI1 and every other UMI1, and between UMI2 and every other UMI2.
    first.row <- all.ind$row <= N
    first.col <- all.ind$col <= N
    keep <- first.row==first.col

    non.flipped.dist <- lev.dist[keep]
    mid <- length(non.flipped.dist)/2L
    to.mid <- seq_len(mid)
    non.flipped.dist <- non.flipped.dist[to.mid] + non.flipped.dist[mid + to.mid]

    non.flipped.row <- true.row[keep][to.mid]
    non.flipped.col <- true.col[keep][to.mid]

    # Pulling out the parts of the matrix corresponding to flipped distances
    # between UMI1 and every other UMI2, and between UMI2 and every other UMI1.
    score.mat <- matrix(NA_integer_, nrow=N, ncol=N)
    low.tri <- !first.row & first.col & true.row > true.col
    score.mat[true.row[low.tri] + (true.col[low.tri]-1L)*N] <- lev.dist[low.tri]
    up.tri <- !first.row & first.col & true.col > true.row
    mat.idx <- true.col[up.tri] + (true.row[up.tri]-1L)*N
    score.mat[mat.idx] <- score.mat[mat.idx] + lev.dist[up.tri]

    # Picking the smallest distance between the flipped and non-flipped distances.
    mat.idx2 <- non.flipped.row + (non.flipped.col-1L)*N
    score.mat[mat.idx2] <- pmin(score.mat[mat.idx2], non.flipped.dist)

    # Making commmunities on those below the distance threshold.
    out <- arrayInd(which(score.mat <= threshold), dim(score.mat))
    g <- make_graph(t(out), n=N)
    return(membership(components(g)))
}
