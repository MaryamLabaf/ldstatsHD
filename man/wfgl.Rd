\name{wfgl}
\alias{wfgl}
\alias{print.wfgl}
\alias{summary.wfgl}

\title{
  weighted fused graphical lasso 
}
\description{

 \code{wfgl} estimates joint partial correlation matrices from two multivariate
  normal distributed datasets using an ADMM based algorithm. Allows for paired data.
}
\usage{

wfgl(D1, D2, lambda1, lambda2, paired = TRUE, automLambdas = TRUE, 
     sigmaEstimate = "CRmad", pairedEst = "Reg-based-sim", maxiter = 30, 
     tol = 1e-05, nsubset = 10000, weights = c(1,1), rho=1, rho.increment = 1, 
     triangleCorrection = TRUE, alphaTri = 0.01, temporalFolders = FALSE, 
     notOnlyLambda2 = TRUE)

}

\arguments{

\item{D1}{
      first population dataset in matrix \eqn{n_1 \times p} form.
}
\item{D2}{
      second population dataset in matrix \eqn{n_2 \times p} form.
}
\item{lambda1}{
      tuning parameter for sparsity in the precision matrices (sequence of lambda1 is allowed).
}
\item{lambda2}{
      tuning parameter for similarity between the precision matrices in the two populations (only one value allowed at a time).
}
\item{paired}{
      if \code{TRUE}, observations in \code{D1} and \code{D2} are assumed to be matched (\eqn{n_1} must be equal to \eqn{n_2}).
}
\item{automLambdas}{
     if \code{TRUE} the lambda's are estimated automatically with \code{lambda1}
     and \code{lambda2} being expected false positive rate levels.
}
\item{sigmaEstimate}{
   robust method used to estimate the variance of estimated partial correlations: name that uniquely identifies  
   \code{"mad"}, \code{"IQR"} or \code{"CRmad"} (default). This measure is used to automatically select the tuning parameter 
   (when \code{automLambdas = TRUE}).
}
\item{pairedEst}{
	type of estimator for the correlation of estimated partial correlation coefficients when \code{"paired = TRUE"}: to select from 
	\code{Reg-based} and \code{Reg-based-sim} (default). This measure is used to weight similarity penalization \code{lambda2} 
	for different pairs of variables.
}
\item{maxiter}{
     maximum number of iterations for the ADMM algorithm.
}
\item{tol}{
     convergence tolerance
}
\item{nsubset}{
     maximum number of estimated partial correlation coefficients (chosen randomly) used to select \code{lambda1}
     and \code{lambda2} automatically (when \code{automLambdas = TRUE}).
}
\item{weights}{
     weights for the two populations to find the inverse covariance matrices.
}
\item{rho}{
     regularization parameter used to compute matrix inverse by eigen value decomposition (default of 1).
}
\item{rho.increment}{
     default of 1.
}
\item{triangleCorrection}{
   if \code{TRUE} the estimated triangle graph structures are tested.
}
\item{alphaTri}{
   significance level for the tested triangle graph structures. 
}
\item{temporalFolders}{
	if \code{TRUE} temporal files are created and eliminated within the procedure. It is used to 
	free R memory space when the dimension is very large (order of thousands).
}
\item{notOnlyLambda2}{
     if \code{FALSE} only lambda2 is found automatically.
}


}
\details{
\code{wfgl} uses a weighted-fused graphical lasso maximum likelihood estimator by solving: 

 \deqn{
 \hat{\Omega}_{WFGL}^{\lambda} = \arg\max\limits_{\Omega_X,\Omega_Y} [\sum_{k=X,Y} \log\det\Omega_k -tr(\Omega_k S_k) - P_{\lambda_1,\lambda_2,V}(\Omega_X,\Omega_Y)],
}
with 
\deqn{
 P_{\lambda_1,\lambda_2, V}(\Omega_X,\Omega_Y) = \lambda_1||\Omega_X||_1 + \lambda_1||\Omega_Y||_1 +\lambda_2\sum_{i,j} v_{ij} |\Omega_{Y_{ij}}-\Omega_{X_{ij}}|,
}
where \eqn{\lambda_1} is the sparsity tuning parameter, \eqn{\lambda_2}  is the similarity tuning parameter, and \eqn{V = [v_{ij}]} 
is a \eqn{p\times p} matrix to weight \eqn{\lambda_2} for each coefficient of the differential precision matrix.
If datasets are independent (\code{paired = "FALSE"}), then it is assumed that \eqn{v_{ij} = 1} for all pairs \eqn{(i,j)}. 
Otherwise (\code{paired = "TRUE"}), weights are estimated in order to account for the dependence structure between datasets in the differential 
network estimation. 

Lambdas can be estimated in each iteration by controlling the expected false positive rate (EFPR) in 
case \code{automLambdas = TRUE}. This transforms the problem of selecting the tuning parameters \eqn{\lambda_1} and \eqn{\lambda_2} to the 
selection of the desired EFPR. In case \code{lambda2} is a single value and \code{lambda1} is a vector with several values, then lambda selection
approaches implemented at \code{\link{lambdaSelection}} can also be used.
  
If \code{triangleCorrection = TRUE}, the weakest edges of estimated triangular motifs are further tested. The reason is that edges that complete 
triangular graph structures suffer an overestimation when applying the ADMM due to using regularized inverse procedures.

}
\value{


  An object of class \code{wfgl} containing the following components:

  \item{path }{adjacency matrices.}
  \item{omega }{precision matrices.}
  \item{triangleCorrection }{determines if triangle structures are tested.}
  \item{weakTriangEdges }{weakest edges in triangle structures which have been tested.}
  \item{weakTriangEdgesPval }{p-values for the weakest edge in triangle structures.}
  \item{diff_value }{convergence control.}
  \item{iters }{number of iterations used.}
  \item{corEst }{dependence structure estimated measure used in the estimation to account for dependent datasets.}
               
}

\references{
Danaher, P., P. Wang, and D. Witten (2014). The joint graphical lasso for inverse covariance estimation across multiple classes. Journal of the Royal Statistical Society: Series B (Statistical Methodology) (2006), 1-20. 

Boyd, S. (2010). Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers. Foundations and Trends in Machine Learning 3(1), 1-122. 
}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}
\seealso{
  \code{\link{plot.wfgl}}  for graphical representation.\cr
  \code{\link{wfrl}}  for weighted fused regression lasso.
}

\examples{
					
# example to use of wfgl
EX2 <- pcorSimulatorJoint(nobs =50, nclusters = 3, nnodesxcluster = c(30, 30,30), 
                          pattern = "pow", diffType = "cluster", dataDepend = "diag", 
                          low.strength = 0.5, sup.strength = 0.9, pdiff = 0.5, nhubs = 5, 
                          degree.hubs = 20,  nOtherEdges = 30, alpha = 2.3, plus = 0, 
                          prob = 0.05, perturb.clust = 0.2, mu = 0, diagCCtype = "dicot", 
                          diagNZ.strength = 0.6, mixProb = 0.5, probSign = 0.7,  
                          exactZeroTh = 0.05)
## not run
#wfgl1 <- wfgl(EX2$D1, EX2$D2, lambda1 = 0.05, lambda2 = 0.1, paired = TRUE, 
#              automLambdas = TRUE, sigmaEstimate = "CRmad", pairedEst = "Reg-based-sim", 
#              maxiter = 30)
#print(wfgl1)

}

