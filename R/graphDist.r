graphDist <- function(A)
{
	  g             <- graph.adjacency(Matrix(A, sparse=T), mode="undirected")
	  g             <- igraph.to.graphNEL(g)
	  DIST          <- as.dist(johnson.all.pairs.sp(g))
	  DIST          <- 1/DIST
	  return(DIST)
}
