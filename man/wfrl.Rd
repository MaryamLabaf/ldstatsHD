\name{wfrl}
\alias{wfrl}
\alias{print.wfrl}
\alias{summary.wfrl}


\title{
  joint estimation of multiple regression coefficient matrices
 }
\description{

	\code{wfrl} estimates jointly two regression coefficient matrices from multivariate
  normal distributed datasets using an ADMM based algorithm.
}
\usage{

wfrl(D1, D2, lambda1, lambda2, automLambdas = TRUE, paired = TRUE, 
	 sigmaEstimate = "CRmad", maxiter=30, tol=1e-05, nsubset = 10000,  
	 rho = 1, rho.increment = 1, notOnlyLambda2 = TRUE)
}

\arguments{
\item{D1}{
      list with the response variables. Two matrices in the list corresponding to the response variables of the two populations.
}
\item{D2}{
      list with the explanatory variables. Two matrices in the list corresponding to the explanatory variables of the two populations.
}
\item{lambda1}{
      tuning parameter for sparsity in the regression coefficients.
}
\item{lambda2}{
      tuning parameter for similarity between the regression coefficients in the two populations.
}
\item{automLambdas}{
     if \code{TRUE} the lambda's are estimated automatically with \code{lambda1}
     and \code{lambda2} being expected false positive rate levels.
}
\item{paired}{
      if \code{TRUE}, observations in \code{D1} and \code{D2} are assumed to be matched (\eqn{n_1} must be equal to \eqn{n_2}).
}
\item{sigmaEstimate}{
   robust method used to estimate the variance of estimated partial correlations: name that uniquely identifies 
   \code{"mad"}, \code{"IQR"} or \code{"CRmad"} (default). This measure is used to automatically select the tuning parameter 
   (when \code{automLambdas = TRUE}).
}
\item{maxiter}{
     maximum number of iterations for the ADMM algorithm.
}
\item{tol}{
     convergence tolerance.
}
\item{nsubset}{
     maximum number of estimated partial correlation coefficients (chosen randomly) used to select \code{lambda1}
     and \code{lambda2} automatically.
}
\item{rho}{
     regularization parameter used to compute matrix inverse by eigen value decomposition (default of 1).
}
\item{rho.increment}{
     default of 1.
}
\item{notOnlyLambda2}{
     if \code{FALSE} only lambda2 is found automatically.
}


}
\details{
\code{wfrl} uses a weighted-fused least squares lasso maximum likelihood estimator by solving: 
 \deqn{
[\hat{\beta}_H,\hat{\beta}_T] =  \arg\min\limits_{\beta_H,\beta_T} \left[ \frac{1}{2n}||Y-\beta_HX||^2_2 + \frac{1}{2n}||Q-\beta_TW||^2_2 +P_{\lambda_1,\lambda_2,V}(\beta)\right]
}
with  
\deqn{
	P_{\lambda_1,\lambda_2,V}(\beta) = \lambda_1||\beta_H||_1 + \lambda_1||\beta_T||_1 + \lambda_2||V \circ (\beta_T-\beta_H)||_1.
}

where \eqn{\lambda_1} is the sparsity tuning parameter, \eqn{\lambda_2}  is the similarity tuning parameter, and \eqn{V = [v_{ij}]} 
is a \eqn{p\times p} matrix to weight \eqn{\lambda_2} for each coefficient of the differential precision matrix.
If datasets are independent (\code{paired = "FALSE"}), then it is assumed that \eqn{v_{ij} = 1} for all pairs \eqn{(i,j)}. 
Otherwise (\code{paired = "TRUE"}), weights are estimated in order to account for the dependence structure between datasets in the differential 
network estimation. An ADMM-type recursive algorithm is used to solve the optimization problem.

See details in  \code{\link{wfgl}} for transforming the selection problem of the tuning parameters \eqn{\lambda_1} and \eqn{\lambda_2}.
}
\value{

An object of class \code{wfrl} containing the following components:

  \item{regCoef }{regression coefficients.}
  \item{path }{ non-zero structure for the regression coefficients.}
  \item{diff_value }{convergence control.}
  \item{iters }{number of iterations used.}

}
\references{
Danaher, P., P. Wang, and D. Witten (2014). The joint graphical lasso for inverse covariance estimation across multiple classes. Journal of the Royal Statistical Society: Series B (Statistical Methodology) (2006), 1-20. 

Boyd, S. (2010). Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers. Foundations and Trends in Machine Learning 3(1), 1-122. 

}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}
\seealso{
   \code{\link{plot.wfrl}}  for graphical representation. \cr
   \code{\link{wfgl}}  for joint partial correlation estimation.
}

\examples{
	
# example to use of wfrl
N	<- 200
EX2 <- pcorSimulatorJoint(nobs = N, nclusters = 3, nnodesxcluster = c(30, 30,30), 
                          pattern = "pow", diffType = "cluster", dataDepend = "diag", 
                          low.strength = 0.5, sup.strength = 0.9, pdiff = 0.5, nhubs = 5, 
                          degree.hubs = 20,  nOtherEdges = 30, alpha = 2.3, plus = 0, 
                          prob = 0.05, perturb.clust = 0.2, mu = 0, diagCCtype = "dicot", 
                          diagNZ.strength = 0.6, mixProb = 0.5, probSign = 0.7,  
                          exactZeroTh = 0.05)
					 
P           <- EX2$P
q           <- 50 
BETA1       <- array(0, dim = c(P, q))
diag(BETA1) <- rep(0.35,q)
BETA2       <- BETA1
diag(BETA2)[c(1:floor(q/2))] <- 0
sigma2      <- 1.3
Q           <- scale(EX2$D1)
W           <- scale(EX2$D2)
X      	    <- Q\%*\%BETA1 + mvrnorm(N,rep(0,q),diag(rep(sigma2,q)))
Y      	    <- W\%*\%BETA2 + mvrnorm(N,rep(0,q),diag(rep(sigma2,q)))
D1     	    <- list(scale(X), scale(Y))
D2     	    <- list(scale(Q), scale(W))
## not run
#wfrl1       <- wfrl(D1, D2, lambda1 = 0.05, lambda2 = 0.05, automLambdas = TRUE, paired = FALSE, 
#                   sigmaEstimate = "CRmad", maxiter = 30, tol = 1e-05, nsubset = 10000, rho = 1, 
#                   rho.increment = 1, notOnlyLambda2 = TRUE)
#print(wfrl1)

}

