\name{plot.wfgl}
\alias{plot.wfgl}

\title{
  Joint partial correlation matrix estimator plot
 }
\description{
 graphical representation of the non-zero joint partial correlation structure.
   
}
\usage{

\method{plot}{wfgl}(x, minn = 0, col = c("blue","red","green"), vertex.size = 3,  
      edgesThickness = FALSE, zoomThick = 10, ...)

}

\arguments{

\item{x}{
      object of class \code{wfgl}.
}
\item{minn}{
      used for visualization purposes in very dense networks. It only plots nodes that have degree larger than \code{minn}.
}
\item{col}{
      vector defining estimated edge colors: common edges (first element), only non-zero coefficients for first population (second element) and
      only non-zero coefficients for second population (third element).
}
\item{vertex.size}{
      \code{\link{plot.igraph}} parameter: vertex sizes.
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


\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}

\seealso{
\code{\link{wfgl}} for joint estimation of multiple partial correlation matrices.
}
\examples{
EX2   <- pcorSimulatorJoint(nobs =50, nclusters = 3, nnodesxcluster = c(30, 30,30), 
                            pattern = "pow", diffType = "cluster", dataDepend = "diag", 
                            diagCCtype = "dicot", diagNZ.strength = 0.6)
## not run
#wfgl1 <- wfgl(EX2$D1, EX2$D2, lambda1 = 0.01, lambda2 = 0.1, paired = TRUE, 
#              automLambdas = TRUE, sigmaEstimate = "CRmad", 
#              pairedEst = "Reg-based-sim", maxiter = 30)
#plot(wfgl1)
#plot(wfgl1, minn = 1, edgesThickness = TRUE)


}


