# UMI filt, sort, group and index.

#UMI_filt=UMI1[id.nfilt]

UMIgroup <- function(UMI_filt, UMI_length = 0, Levensh_threshold = 2)
{
    # Sorting UMI and grouping by levenshtein distance.
    UMI_sort <- DNAStringSet(names(sort(table(UMI_filt), decreasing=TRUE)))
    UMI_levenshtein.out <- .UMIlevenshtein(UMI_sort, levenshtein_max = UMI_length)

    sortedUMI <- vector("list", length(UMI_levenshtein.out))
    for (i in 1:length(UMI_levenshtein.out)){
        if (length((which(UMI_levenshtein.out[[i]]<=Levensh_threshold&sort(table(UMI_filt), decreasing=TRUE)[[i]]>=2*sort(table(UMI_filt), decreasing=TRUE)-1)))!=0){
            sortedUMI[[i]] <- rbind(i, which(UMI_levenshtein.out[[i]]<=Levensh_threshold&sort(table(UMI_filt), decreasing=TRUE)[[i]]>=2*sort(table(UMI_filt), decreasing=TRUE)-1))
        }else{
            sortedUMI
        }
    }
    
    # Connected UMIs are shown in the matrix and used to make a graph.
    sortedUMI.bind <- do.call(cbind, sortedUMI)
    if(!is.null(sortedUMI.bind)){
        UMI_mat_prep <- matrix(sortedUMI.bind, nc = 2, byrow = TRUE)
        UMI_mat_prep_rev <- cbind(UMI_mat_prep[,2], UMI_mat_prep[,1])

        for(i in 1:nrow(UMI_mat_prep)){
            for(a in 1:nrow(UMI_mat_prep)){ 
                if (sum(UMI_mat_prep[i,]==UMI_mat_prep_rev[a,]) > 1){
                    UMI_mat_prep[a,] <- c(0,0)
                }
            }
        }
        UMI_mat_prep <- UMI_mat_prep[-c(which(UMI_mat_prep[,1]==0)),]
        UMI_graph <- make_graph(as.vector(t(UMI_mat_prep)), n = length(sortedUMI), directed = FALSE)

        # Indexing the input UMIs to be able to subset the corresponding reads.
        UMI_cluster <- split(UMI_sort, membership(components(UMI_graph)))

        UMI_index <- vector("list", length(UMI_cluster))
        for (i in 1:length(UMI_cluster)){
            len <- 1
            for (a in 1:length(UMI_cluster[[i]])){
                id <- which(UMI_filt==UMI_cluster[[i]][a])
                for (b in 1:length(id)){
                    UMI_index[[i]][len] <- id[b] 
                    len <- len+1
                }
            }
        }
    
    }else{
        UMI_index <- NA
        
    }
    
    return(UMI.id = UMI_index)
}

# Calculates levenshtein distance between UMIs.
.UMIlevenshtein <- function(UMI_sort, levenshtein_max = 6)
{
    UMI_levenshtein <- vector("list", length(UMI_sort))
    
    for (i in 1:length(UMI_sort)){
        UMI_levenshtein[[i]] <- stringdist(UMI_sort[i], UMI_sort)
        UMI_levenshtein[[i]][i] <- levenshtein_max
    } 
    
    return(UMI_levenshtein)
}