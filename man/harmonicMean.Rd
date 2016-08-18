\name{harmMeanDist}
\alias{harmMeanDist}


\title{
  Harmonic mean of network distances
 }
\description{

	Finds the harmonic mean of the geodesic distances between nodes in a graph.
	
}
\usage{
harmMeanDist(A, nodesDegree = NULL)

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
 The geodesic distance \eqn{d_{ij}} between two nodes \eqn{i} and \eqn{j} of a graph is 
 given by the shortest number of edges so that we can go from one of the nodes
 to the other. Then, the harmonic mean of these distances is given by 
\deqn{
\bar{d} = [p(p-1)/(2\sum_{i<j} d_{ij})]^{-1}.
}

}
\value{          
 The harmonic mean value.
}
\references{
Costa, L. and F. Rodrigues (2007). Characterization of complex networks: A survey of measurements. Advances in Physics 56(1), 167-242. 
}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}
\seealso{
 \code{\link{lambdaSelection}} for lambda selection approaches.\cr
}

\examples{

# example to use of harmonic mean function
EX1         <- pcorSimulator(nobs = 50, nclusters = 2, nnodesxcluster = c(40,30), 
                             pattern = "powerLaw")
y           <- EX1$y
out3        <- huge(y, method = "mb", lambda = 0.4)
PATH        <- out3$path[[1]]
hm          <- harmMeanDist(PATH)
  
 
}


