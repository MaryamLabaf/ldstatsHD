\name{graphDist}
\alias{graphDist}


\title{
  Inverse of the geodesic distance matrix
 }
\description{

	\code{graphDist} computes the dissimilarity matrix (inverse of geodesic distance) of a graph structure.
	
}
\usage{

graphDist(A)
}
\arguments{

\item{A}{
      \code{\link{matrix}} or \code{\link{Matrix}} object with adjacency matrix of a graph.
}

}
\details{
 The geodesic distance between two nodes \eqn{i} and \eqn{j} of a graph is 
 given by the shortest number of edges so that we can go from one of the nodes
 to the other.

}
\value{          
	a dissimilarity lower triangular matrix with the inverse of the geodesic distance matrix.

}
\references{
Costa, L. and F. Rodrigues (2007). Characterization of complex networks: A survey of measurements. Advances in Physics 56(1), 167-242. 
}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}
\seealso{
 \code{\link{graphCorr}} for graph correlation matrix.\cr
}

\examples{
# example to use of graphDist function
EX1         <- pcorSimulator(nobs = 50, nclusters = 2, nnodesxcluster = c(40,30), 
                             pattern = "powerLaw")
y           <- EX1$y
out3        <- huge(y, method = "mb", lambda = 0.4)
gd          <- graphDist(out3$path[[1]])
   
 
}


