harmMeanDist <- function(A, nodesDegree = NULL)
{
    if(is.null(nodesDegree))   nodesDegree    <- degrees(A)
    P		  <- dim(A)[1]
    minDegree <- 2
    
    if (any(nodesDegree >= minDegree))
    {
      NODESINGRAPH  <- which(nodesDegree > 0)
      g             <- graph.adjacency(Matrix(A[NODESINGRAPH,NODESINGRAPH], sparse=T),mode="undirected")
      g             <- igraph.to.graphNEL(g)
      DIST          <- as.dist(johnson.all.pairs.sp(g))
      WHICH.DIST    <- P * (P-1) / (sum(1/DIST[is.finite(DIST)]) * 2)
    }
    else
     WHICH.DIST <- 0
    return(WHICH.DIST)
}


