\name{graphCorr}
\alias{graphCorr}


\title{
  One minus graph correlation matrix
 }
\description{

	\code{graphCorr} computes the dissimilarity matrix (one minus the correlation 
  matrix) of a graph structure.
	
}
\usage{

graphCorr(A, nodesDegree = NULL)
}
\arguments{

\item{A}{
      \code{\link{matrix}} or \code{\link{Matrix}} object with adjacency matrix of a graph.
}
\item{nodesDegree}{
      vector with nodes degree (in case it is been previously calculated).
}

}
\details{
The similarity matrix of a graph is given by
\deqn{
\sigma_{ij} = \frac{\eta_{ij}}{\sqrt{k_ik_j}},
}
where \eqn{\eta_{ij}} is defined by the number of common neighbors of nodes \eqn{i} and \eqn{j}
and \eqn{k_i} is the degree of the node \eqn{i}. The dissimilarity matrix is given by
\deqn{
\delta_{ij} = 1- \sigma_{ij}.
}
 

}
\value{          
	a dissimilarity lower triangular matrix with one minus the correlation of the graph nodes.

}
\references{
Costa, L. and F. Rodrigues (2007). Characterization of complex networks: A survey of measurements. Advances in Physics 56(1), 167-242. 
}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}
\seealso{
 \code{\link{graphDist}} for inverse of the geodesic distance matrix.\cr
}

\examples{

# example to use of graphCorr function
EX1         <- pcorSimulator(nobs = 50, nclusters = 2, nnodesxcluster = c(40,30), 
                             pattern = "powerLaw")
y           <- EX1$y
out3        <- huge(y, method = "mb", lambda = 0.4)
gc          <- graphCorr(out3$path[[1]])
  
 
}


