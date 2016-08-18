\name{plot.pcorSim}
\alias{plot.pcorSim}

\title{
  Partial correlation matrix simulator plot
 }
\description{
 graphical representation of the non-zero partial correlation matrix structure
   
}
\usage{

\method{plot}{pcorSim}(x, vertex.size = 3, vertex.label = NA, hubsCol = TRUE, ...)

}

\arguments{

\item{x}{
      object of class \code{pcorSim}.
}
\item{vertex.size}{
      \code{\link{plot.igraph}} parameter: vertex sizes.
}
\item{vertex.label}{
      \code{\link{plot.igraph}} parameter: vertex label.
}
\item{hubsCol}{
      if \code{TRUE} hub nodes are highlighted in a different color.
}
\item{\dots}{
    arguments passed to or from other methods to the low level.
}

}

\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}

\seealso{
\code{\link{pcorSimulator}} for partial correlation matrix generation.\cr
\code{\link{pcorSimulatorJoint}} for joint partial correlation matrix generation.
}
\examples{

EX1 <- pcorSimulator(nobs = 50, nclusters = 3, nnodesxcluster = c(100,30,50), 
                     pattern="powerLaw", plus = 0)
plot(EX1)
                
EX2 <- pcorSimulator(nobs = 25, nclusters = 2, nnodesxcluster = c(60,40), 
                     pattern = "powerLaw", plus = 1)
plot(EX2)
 
}

