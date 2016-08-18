\name{plot.eqCorTestByRows}
\alias{plot.eqCorTestByRows}

\title{
   plot for equality of two correlation matrices by rows test  
 }
\description{
 graphical representation for the equality of two correlation matrices test by rows: confidence intervals of the test statistics.   
}
\usage{

\method{plot}{eqCorTestByRows}(x, mains = c("AS CI", "max CI"), xlabs = c("",""), ylabs = c("",""), 
      pch = "-", ownCols = TRUE, ...)
}

\arguments{

\item{x}{
      object of class \code{eqCorTestByRows}.
}
\item{mains}{
      vector of size two with main of plots (for average of squares test and extreme value test).
}
\item{xlabs}{
      vector of size two with xlabs of plots (for average of squares test and extreme value test).
}
\item{ylabs}{
      vector of size two with ylabs of plots (for average of squares test and extreme value test).
}
\item{pch}{
  	  \code{pch} given to identify confidence interval limits.
}
\item{ownCols}{
  	  if \code{ownCols = TRUE} green and black colors are generated with green lines identifying significant variables.
}
\item{\dots}{
    arguments passed to or from other methods to the low level.
}

}

\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}

\seealso{
\code{\link{eqCorTestByRows}} for equality of two correlation matrices by rows test. 
}
\examples{

EX2 <- pcorSimulatorJoint(nobs = 200, nclusters = 3, nnodesxcluster = c(60,40,50), 
                         pattern = "pow", diffType = "cluster", dataDepend = "diag", 
                         pdiff=0.5)
			               
#### eq corr by rows
## not run
#eqCorEX2 <- eqCorTestByRows(EX2$D1, EX2$D2, testStatistic = c("AS", "max"), 
#                            nite = 200, paired = TRUE, exact = TRUE, 
#                            subMatComp = FALSE, iniP = 1, finP = 40, 
#                            conf.level = 0.95)
#plot(eqCorEX2)         					

 
}

