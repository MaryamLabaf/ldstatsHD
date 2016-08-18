\name{plot.pcorSimJoint}
\alias{plot.pcorSimJoint}

\title{
  Joint partial correlation matrix simulator plot
 }
\description{
 graphical representation of two non-zero partial correlation structures
   
}
\usage{

\method{plot}{pcorSimJoint}(x, minn = 0, col = c("blue","red","green"), vertex.size=3, 
      edgesThickness = FALSE, ...)

}

\arguments{

\item{x}{
      object of class \code{pcorSimJoint}.
}
\item{minn}{
      used for visualization purposes in very dense networks. It only plots nodes that have degree larger than \code{minn}.
}
\item{col}{
      vector defining edge colors: common edges (first element), only non-zero coefficients for first population (second element) and
      only non-zero coefficients for second population (third element).
}
\item{vertex.size}{
      \code{\link{plot.igraph}} parameter: vertex sizes.
}
\item{edgesThickness}{
  	  if \code{TRUE}, an edge thickness is proportional to the magnitude of its underlying estimated partial correlation coefficient.
}
\item{\dots}{
    arguments passed to or from other methods to the low level.
}

}

\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}

\seealso{
\code{\link{pcorSimulatorJoint}} for joint partial correlation matrix generation.
}
\examples{
EX1 <- pcorSimulatorJoint(nobs = 50, nclusters = 2, nnodesxcluster = c(30, 40), 
                         pattern = "pow", diffType = "cluster", dataDepend = "ind", 
                         pdiff = 0.2, diagCCtype = "dicot", diagNZ.strength = .5)
plot(EX1, edgesThickness = TRUE)

EX2 <- pcorSimulatorJoint(nobs = 50, nclusters = 3, nnodesxcluster = c(30, 40,60), 
                          pattern = "pow", diffType = "cluster", dataDepend = "diag", 
                          pdiff = 0.4, diagCCtype = "beta")
plot(EX2)

 
}


