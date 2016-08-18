\name{agnesCoef}
\alias{agnesCoef}


\title{
  AGNES coefficient 
 }
\description{

	computes the AGNES coefficient of a graph structure.
	
}
\usage{
	
agnesCoef(A)

}
\arguments{

\item{A}{
      \code{\link{matrix}} or \code{\link{Matrix}} object with adjacency matrix of a graph.
}

}
\details{
The input of the AGNES hierarchical algorithm is the dissimilarity matrix of the graph structure \code{A} 
computed by \code{\link{graphCorr}}. Then the R function \code{\link{agnes}} is used to calculate the coefficient.
}
\value{          
AGNES coefficient magnitude for the given graph structure. 

}
\references{
Rousseeuw, P., A. Struyf, and M. Hubert (2013). cluster: Cluster Analysis Basics and Extensions.
}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}
\seealso{
 \code{\link{lambdaSelection}} for lambda selection approaches.\cr
}

\examples{

EX1        <- pcorSimulator(nobs = 50, nclusters = 2, nnodesxcluster = c(40,30), 
                            pattern = "powerLaw")
y           <- EX1$y
out3        <- huge(y, method = "mb", lambda = 0.4)
PATH        <- out3$path[[1]]
hm          <- agnesCoef(PATH)

 
}


