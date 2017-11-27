# Function for both UMIs.
# id.filt needs to be specified.
# umiGroup adds the levenshtein distances of both UMI's up and then proceeds as in the initial UMIgroup function.

umiGroup <- function(UMI1, UMI2 = NULL, groups=NULL, threshold=2) {
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
        UMI.id <- .umiClust(UMI1=UMI1[current], UMI2=UMI2[current], threshold=threshold)
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
    
    to.link <- which(lev.dist <= threshold)
    if (length(to.link)==0L) { 
        return(seq_along(UMI1))
    }

    # Getting indices without expanding the dist structure.
    N <- length(UMI1)
    col.ends <- cumsum((N-1L):1L)
    current.col <- findInterval(to.link, col.ends, left.open=TRUE) + 1L
    row.starts <- seq_len(N-1L)
    current.row <- row.starts[current.col] + to.link - c(0L, col.ends)[current.col]  

    # Identifying the components of the graph.
    g <- make_graph(rbind(current.col, current.row), n=N)
    return(membership(components(g)))
}

