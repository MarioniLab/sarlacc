UMIgraph <- function(UMImat)

#The assigned connections between the UMI's are used as edges to construct a graph.

{
    UMI_graph <- graph.adjacency(UMImat)
    plot(UMI_graph)
    
    return(UMI_graph)
}
