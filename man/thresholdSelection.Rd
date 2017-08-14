\name{thesholdSelection}
\alias{thesholdSelection}


\title{
   Threshold selection for equality of correlation matrices test
 }
\description{

This function finds the threshold that maximizes an estimate of the power of the exceedances-based test.

}
\usage{

thesholdSelection(D1, D2, useq, deltaA = 3, deltaB = 10, nite = 500, 
                  excAdj = FALSE, alpha = 0.05, paired = TRUE)
                
}

\arguments{

\item{D1}{
      first population dataset in matrix \eqn{n\times p} form.
}
\item{D2}{
      second population dataset in matrix \eqn{n\times p} form.
}
\item{useq}{
      sequence of threshold levels to be used.
}
\item{deltaA}{
	 shape hyperparameter for gamma prior distribution of Fisher transform correlation differences in absolute value.
}
\item{deltaB}{
	 scale hyperparameter for gamma prior distribution of Fisher transform correlation differences in absolute value.
}
\item{nite}{
  	 number of generated samples.
}
\item{excAdj}{
      weight for the exceedances test. If \code{excAdj = FALSE} the test statistic is given by the squared exceedances. In contrast, If \code{excAdj = TRUE} 
      the test statistic is given by the squared of the exceedances minus the threshold.
}
\item{alpha}{
      null hypothesis rejection level.
}
\item{paired}{
      if \code{TRUE}, observations in \code{D1} and \code{D2} are assumed to be matched (\eqn{n_1} must be equal to \eqn{n_2}).
}

\item{\dots}{
    arguments passed to or from other methods to the low level.
}

}
\details{
 Details are given in references.
  }
\value{
	Optimal threshold
}
\references{
To come
}

\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina, Claus Mayer and Ioannis Papastathopoulos. 
}

\examples{
#### data
EX2 <- pcorSimulatorJoint(nobs = 50, nclusters = 3, nnodesxcluster = c(40,40,40), 
                          pattern = "pow", diffType = "cluster", dataDepend = "diag", 
                          pdiff=0.5)			               
						 
#### best threshold
useq <- seq(0,qnorm(1-0.01/2),length.out=150)
thesholdSelection(EX2$D1,EX2$D2, useq)

}
