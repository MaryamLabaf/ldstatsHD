\name{eqCorrMatTest}
\alias{eqCorrMatTest}
\alias{print.eqCorrMatTest}
\alias{summary.eqCorrMatTest}


\title{
   equality of two correlation matrices test
 }
\description{

	tests the equality of two correlation matrices coming from two independent (or paired) datasets, 
	that can possibly be high dimensional. 
}
\usage{

eqCorrMatTest(D1, D2 = NULL, testStatistic = c("AS", "max", "exc"), 
              testNullDist = c("asyIndep","asyDep", "np"), nite= 500, 
              paired = FALSE, threshold = 2.3, excAdj = TRUE, exact = FALSE,
              conf.level = 0.95, saddlePoint = FALSE, MINint =2, MAXint=100, ...)
                
}

\arguments{

\item{D1}{
      first population dataset in matrix \eqn{n_1\times p} form.
}
\item{D2}{
      second population dataset in matrix \eqn{n_2\times p} form. If \code{D2 = NULL}, an hypothesis testing for non-identity 
      correlation matrix is performed instead.
}
\item{testStatistic}{
      test statistic used for the hypothesis testing: name that uniquely identifies  \code{"AS"} for average of squares
      based test, \code{"max"} for an extreme values test and \code{"exc"} for
      an exceedances-based test. If \code{"all"}, all tests are performed.
}
\item{testNullDist}{
      Null distribution approximation. Either assuming independence between sample coefficients (\code{"asyIndep"}), accounting for dependence in a 
      parametric distribution (\code{"asyDep"}), or accounting for dependence by a non-parametric distribution (\code{"np"}).
}
\item{nite}{
      number of iterations used to generate the permuted samples (only if 
      \code{testNullDist = "asyDep"} or \code{testNullDist = "np"}). It is also used to approximate a null distribution by Monte Carlo when 
      \code{testNullDist = "asyIndep"} and  \code{testStatistic = "exc"}. 
}
\item{paired}{
      if \code{TRUE}, observations in \code{D1} and \code{D2} are assumed to be matched (\eqn{n_1} must be equal to \eqn{n_2}).
}
\item{threshold}{
      exceedance threshold used if \code{testStatistic = "exc"}.
}
\item{excAdj}{
      weight for the exceedances test. If \code{excAdj = FALSE} the test statistic is given by the squared exceedances. In contrast, If \code{excAdj = TRUE} 
      the test statistic is given by the squared of the exceedances minus the threshold.
}
\item{exact}{
      permuted samples method: if \code{TRUE} it forces to have the exact same number of observations 
      in the two conditions in the samples exchanging process. If \code{FALSE}, permutations are made exchanging 
      matched observations from the two datasets randomly with probability equal to 0.5.
}
\item{conf.level}{
	confidence level of the interval.
}
\item{saddlePoint}{
	saddle-point approximation is used to find the p-value for exceedances-based test under asymptotic independence null distribution (only for \code{excAdj = FALSE}). 
}
\item{MINint}{
	minimum expected number of exceedances for saddle point approximation. 
}
\item{MAXint}{
	maximum expected number of exceedances for saddle point approximation. 
}
\item{\dots}{
    arguments passed to or from other methods to the low level.
}

}
\details{
  The extreme value test is the most powerful test against very sparse alternatives whereas
  sum of squares test is the most powerful when the true differential correlation 
  matrix is dense. Otherwise, for a reasonable selection of the 
  exceedance threshold, the exceedances test overperforms the power of the other two tests.
  
  Paired structures can be used in this function. For instance, if \code{paired = TRUE} then the correlation
  between sample correlation coefficients in the two matrices is estimated to adjust the test statistic.
  
  Asymptotic independence tests are fast since they do not compute permuted samples and can be used, even for paired data, under
  weakly dependent assumptions (very sparse correlation matrices) when sample sizes are large. If these assumptions do not hold, 
  wrong representations of the p-values under \eqn{H_0} could be found, in which case, permuted based null distributions should be used instead.
  
  Testing if a correlation matrix is the identity matrix can also be performed when \code{D2 = NULL}. Note that the same type of test statistics
  and null distributions are available in this setting. Nevertheless, Monte Carlo simulations are used instead of permuted samples. Here null distributions are approximated by computing sample correlation matrices of generated data following a \eqn{N_p(0,I)} under the assumption of normality.
 }
\value{
	An object of class \code{eqCorrMatTest} containing the test statistics, p-values and confidence intervals for the selected tests.
}
\references{
To come
}

\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina, Claus Mayer and Ioannis Papastathopoulos. 
}
\seealso{
\code{\link{eqCorTestByRows}} for testing linear independence and equality of sames rows in two correlation matrices.\cr
}

\examples{
#### data
EX2 <- pcorSimulatorJoint(nobs = 50, nclusters = 3, nnodesxcluster = c(40,40,40), 
                          pattern = "pow", diffType = "cluster", dataDepend = "diag", 
                          pdiff=0.5)			               
						 
#### eq. corr. mat. test
## not run
#test1 <- eqCorrMatTest(EX2$D1, EX2$D2, testStatistic = c("AS", "exc", "max"), 
#                       paired = TRUE, nite = 400)
#print(test1)

## not run
#test2 <- eqCorrMatTest(EX2$D1, NULL, testStatistic = c("AS", "exc", "max"), 
#                       paired = TRUE, nite = 400)
#print(test2)

}
