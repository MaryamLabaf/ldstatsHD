\name{plot.wfrl}
\alias{plot.wfrl}

\title{
  Joint regression coefficient matrix estimator plot
 }
\description{
 graphical representation of the non-zero joint regression coefficients structure
   
}
\usage{

\method{plot}{wfrl}(x, minn = 0, col = c("blue","red","green"), vertex.size = 2, 
      vertex.color = c("red", "blue"), edgesThickness = FALSE, 
      zoomThick = 10, ...)

}

\arguments{

\item{x}{
      object of class \code{wfrl}.
}
\item{minn}{
      used visualization purposes in very dense networks. It only plots nodes that have degree larger than \code{minn}.
}
\item{col}{
      vector defining estimated edge colors: common edges (first element), only non-zero coefficients for first population (second element) and
      only non-zero coefficients for second population (third element).
}
\item{vertex.size}{
      \code{\link{plot.igraph}} parameter: vertex sizes.
}
\item{vertex.color}{
      vector defining the vertex colors for directed graph: first element describes the color of explanatory variables and 
      second element describes the color for response variables.
}
\item{edgesThickness}{
  	  if \code{TRUE}, an edge thickness is proportional to the magnitude of its underlying estimated partial correlation coefficient.
}
\item{zoomThick}{
	it increases the thickness of all edges by \code{zoomThick} times (used for visualization purposes).
}
\item{\dots}{
    arguments passed to or from other methods to the low level.
}

}
\details{
	It produces a directed graph structure that connects explanatory variables to response variables.
}

\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}

\seealso{
\code{\link{wfrl}} for joint estimation of regression coefficients.\cr
}
\examples{
N	<- 200
EX2 <- pcorSimulatorJoint(nobs = N, nclusters = 3, nnodesxcluster = c(60,40,50), 
                          pattern = "pow", diffType = "cluster", dataDepend = "diag", 
                          low.strength = 0.5, sup.strength = 0.9, pdiff = 0.5, nhubs = 5, 
                          degree.hubs = 20,  nOtherEdges = 30, alpha = 2.3, plus = 0, 
                          prob = 0.05, perturb.clust = 0.2, mu = 0, diagCCtype = "dicot", 
                          diagNZ.strength = 0.6, mixProb = 0.5, probSign = 0.7,  
                          exactZeroTh = 0.05)
					 
P           <- EX2$P
q           <- 50 
BETA1       <- array(0,dim=c(P,q))
diag(BETA1) <- rep(0.35,q)
BETA2       <- BETA1
diag(BETA2)[c(1:floor(q/2))]<-0
sigma2 	<- 1.3
Q       <- scale(EX2$D1)
W      	<- scale(EX2$D2)
X      	<- Q\%*\%BETA1 + mvrnorm(N,rep(0,q),diag(rep(sigma2,q)))
Y      	<- W\%*\%BETA2 + mvrnorm(N,rep(0,q),diag(rep(sigma2,q)))
D1     	<- list(scale(X),scale(Y))
D2     	<- list(scale(Q),scale(W))
## not run
#wfrl1   <- wfrl(D1, D2, lambda1 = 0.01, lambda2 = 0.05, automLambdas = TRUE, paired = FALSE, 
#               sigmaEstimate = "CRmad", maxiter=30, tol=1e-05)
#plot(wfrl1)
#plot(wfrl1, minn = 1, edgesThickness = TRUE)

}


