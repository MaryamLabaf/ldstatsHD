\name{eqCorTestByRows}
\alias{eqCorTestByRows}
\alias{print.eqCorTestByRows}
\alias{summary.eqCorTestByRows}

\title{
  Correlation matrices test by rows
 }
\description{
 Tests whether the \eqn{g}th row of a correlation matrix is either non-zero or different to the same row of another correlation matrix. Allows for paired data. 
}

\usage{
eqCorTestByRows(D1, D2 = NULL, testStatistic = c("AS", "max"), nite = 200, 
                paired = FALSE, exact = TRUE, subMatComp = FALSE, iniP = 1, 
                finP = NULL, conf.level = 0.95)

}

\arguments{

\item{D1}{
      first population dataset in matrix \eqn{n_1\times p} form.
}
\item{D2}{
      second population dataset in matrix \eqn{n_2\times p} form. If \code{D2 = NULL} non-zero correlation rows test is performed instead.
}
\item{testStatistic}{
      test statistic used for the hypothesis testing: name that uniquely identifies  \code{"AS"} for average of squares
      based test and \code{"max"} for an extreme value test.
}
\item{nite}{
      number of iterations used to generate the permuted samples.
}
\item{paired}{
      if \code{TRUE}, observations in \code{D1} and \code{D2} are assumed to be matched (\eqn{n_1} must be equal to \eqn{n_2}).
}
\item{exact}{
      permuted samples method: if \code{TRUE} it forces to have the exact same number of observations 
      in the two conditions in the samples exchanging process. If \code{FALSE}, permutations are made exchanging 
      matched observations from the two datasets randomly with probability equal to 0.5.
}
\item{subMatComp}{
	 used to reduce computational time when using the test in very high dimensional data. If \code{TRUE} correlation sub-matrices are 
	 used, if \code{FALSE}, whole correlation matrices are computed.
}
\item{iniP}{
      only for \code{subMatComp = TRUE} and \code{D2 != NULL}. First row to be tested.
}

\item{finP}{
      only for \code{subMatComp = TRUE} and \code{D2 != NULL}. Last row to be tested.
}
\item{conf.level}{
	confidence level of the interval.
}

}
\details{
	This test uses a sum of squares based test statistic as given by the adjusted squared correlation \code{\link{cor2mean.adj}} 
	as well as an extreme value based test statistic as given by \code{\link{max}}. 
		
	Null distributions are approximated differently when  testing equality of two correlation rows and testing if correlation rows are equal to zero. 
	In the first case, permuted samples are used to construct the confidence interval (see details in \code{\link{eqCorrMatTest}}).  
	In the latter, they are found using Monte Carlo samples. For instance, \eqn{n} iid observations from a normal distribution \eqn{N(0,1)} are generated. 
	Then, the adjusted square (or absolute maximum) correlations between these montecarlo samples and the original data \eqn{D1} are found. 
}

\value{
  An object of class \code{eqCorTestByRows} containing the following components:
  \item{AStest }{average of squares test statistics.}
  \item{pvalAS }{average of squares test p-values.}
  \item{ciAS }{average of of squares test statistic confidence interval.}
  \item{Maxtest }{extreme value test statistics.}
  \item{pvalMax }{extreme value test p-values.}
  \item{ciMax }{extreme value test statistic confidence interval.}
}
\references{
to come.
}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}
\seealso{
\code{\link{plot.eqCorTestByRows}} for graphical representation.\cr
\code{\link{eqCorrMatTest}} for testing equality of two correlation matrices.

}

\examples{
#### data
EX2 <- pcorSimulatorJoint(nobs = 200, nclusters = 3, nnodesxcluster = c(60,40,50), 
                          pattern = "pow", diffType = "cluster", dataDepend = "diag", 
                          pdiff = 0.5)
			               
#### eq corr by rows
## not run
#test1 <- eqCorTestByRows(EX2$D1, EX2$D2, testStatistic = c("AS", "max"), 
#                         nite = 200, paired = TRUE, exact = TRUE, subMatComp = FALSE, 
#                         iniP = 1, finP = 40, conf.level = 0.95)
#print(test1)

#### zero corr by rows
#test2 <- eqCorTestByRows(EX2$D1, testStatistic = c("AS", "max"), nite = 1000, 
#                         conf.level = 0.95)
#print(test2)

}


