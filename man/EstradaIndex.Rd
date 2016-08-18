\name{estradaIndex}
\alias{estradaIndex}


\title{
  Estrada Index of a graph structure
 }
\description{

	Computes the Estrada Index given the adjacency matrix of a graph structure.
	
}
\usage{
estradaIndex(A)
}
\arguments{

\item{A}{
      \code{\link{matrix}} or \code{\link{Matrix}} object with adjacency matrix of a graph.
}

}
\details{
The Estrada Index is calculated by
\deqn{
EE(\lambda)  = \sum_{j=1}^p \exp(\gamma_j(\lambda)),
}
where \eqn{\gamma_1(\lambda), \ldots,\gamma_p(\lambda)} are the eigenvalues of \eqn{A_G^\lambda}. 
 

}
\value{          
	Estrada index coefficient
}
\references{
Estrada, E. (2011). The structure of complex networks. New York: OXFORD University press. 
}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}
\seealso{
 \code{\link{lambdaSelection}} for lambda selection approaches.\cr
}

\examples{

EX1         <- pcorSimulator(nobs = 50, nclusters=2, nnodesxcluster=c(40,30), 
                             pattern="powerLaw")
y           <- EX1$y
out3        <- huge(y, method = "mb", lambda = 0.4)
PATH        <- out3$path[[1]]
hm          <- estradaIndex(PATH)
  
 
}

