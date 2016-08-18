\name{pcorSimulatorJoint}
\alias{pcorSimulatorJoint}
\alias{print.pcorSimulatorJoint}
\alias{summary.pcorSimulatorJoint}

\title{
  joint partial correlation matrices simulator 
 }
\description{

	\code{pcorSimulatorJointPaired} creates two similar positive definite precision matrices with three possible 
   graph structures: hubs-based, power-law and random. Moreover, it allows for three types of 
   differential graph structures: random differences, clustered differences or a mixture of the two. Then, it 
   generates (dependent) datasets from a multivariate normal distribution defined by the inverse of such precision matrices.
   
}
\usage{

pcorSimulatorJoint(nobs, nclusters, nnodesxcluster, pattern = "hubs", 
                   diffType = "cluster", dataDepend = "ind", low.strength = 0.5, 
                   sup.strength = 0.9, pdiff = 0, nhubs = 5, degree.hubs = 20,  
                   nOtherEdges = 30, alpha = 2.3, plus = 0, prob = 0.05, 
                   perturb.clust = 0, mu = 0, diagCCtype = "dicot", 
                   diagNZ.strength = .5, mixProb = 0.5, probSign = 0.5,  
                   exactZeroTh = 0.05, seed = sample(10000,nclusters+2))
                  
}

\arguments{

\item{nobs}{
      number of observations.
}
\item{nclusters}{
      number of clusters or blocks of variables.
}
\item{nnodesxcluster}{
      number of nodes/variables per cluster.
}
\item{pattern}{
      graph structure pattern: name that uniquely identifies  \code{"hubs"}, \code{"power"} and \code{"random"}.
}
\item{diffType}{
     pattern in differential edges: name that uniquely identifies  \code{"random"}, \code{"cluster"}  or \code{"mixed"}.
}
\item{dataDepend}{
 model used to describe the dependent structure for the data: name that uniquely identifies  \code{"ind"} (no dependence), \code{"diagOmega"}, \code{"mult"} or \code{"add"}.
}
\item{low.strength}{
      minimum magnitude for nonzero partial correlation elements before regularization.
}
\item{sup.strength}{
      maximum magnitude for nonzero partial correlation elements before regularization.
}
\item{pdiff}{
     proportion of differential edges from the total number edges in each graph.
}
\item{nhubs}{
      number of hubs per cluster (if \code{pattern = "hubs"}).
}
\item{degree.hubs}{
      degree of hubs (if \code{pattern = "hubs"}).
}
\item{nOtherEdges}{
      number of edges for non-hub nodes (if \code{pattern = "hubs"}).
}
\item{alpha}{
      positive coefficient for the Riemman function in power-law distributions.
}
\item{plus}{
      power-law distribution added complexity (zero by default).
}
\item{prob}{
      probability of edge existence for random networks (if \code{pattern="random"}).
}
\item{perturb.clust}{
     proportion of the total number of edges that are connecting two different clusters.
}
\item{mu}{
     expected values vector to generate data (zero by default).
}
\item{diagCCtype}{
 way to generate diagonal values of either cross partial correlation matrix (if \code{dataDepend = "diagOmega"}) or cross correlation matrix (if \code{dataDepend = "mult"} or \code{dataDepend = "add"}): name that uniquely identifies  \code{"dicot"} or \code{"beta13"} (see details).
}
\item{diagNZ.strength}{
 magnitude for the non-zero elements in the diagonal of the cross (partial) correlation when \code{diagCCtype = "dicot"}.
}
\item{mixProb}{
     proportion of random differential connections if \code{diffType = "mixed"}. The remaining
     connections are given by a cluster type.
}
\item{probSign}{
     probability of positive sign for non-zero partial correlation coefficients. Thus, negative signs 
     are obtained with probability \code{1-probSign}. 
}
\item{exactZeroTh}{
     partial correlation coefficients smaller than exactZeroTh are considered exact zeros.
}
\item{seed}{
      vector with seeds for each cluster.
}
}
\details{
First, \code{\link{pcorSimulator}} is used to create a common precision matrix among the two populations. 
Then, differential edges are added based on the next two patterns: Cluster - a graph cluster is zero in one condition and non-zero in the other   
condition; Random - differential connections are given randomly in the graph. 
  
Paired structure is defined by arguments \code{dataDepend} and \code{diagCCtype}. 
Additive (\code{dataDepend = "add"}) and multiplicative (\code{dataDepend = "mult"}) models are used on the cross-covariance matrix
such that \eqn{\Sigma_{XY} = \Delta \Sigma_X \Delta^t}, with diagonal matrix \eqn{\Delta}, \eqn{0\leq\Delta_{ii}<1} and \eqn{\Sigma_{XY} = \Delta\Sigma_X^{1/2}\Sigma_Y^{1/2} \Delta^t} respectively where diagonal coefficients in \eqn{\Delta} are defined by \code{diagCCtype}. 
A simplification is also considered by assuming that variables in one data set are only
conditionally dependent to the same variables of the other data set, hence assuming
a diagonal structure in the cross joint partial correlation matrix that can also be defined by \eqn{\Delta}. For the three models,
In case \code{diagCCtype = "dicot"} the diagonal elements in \eqn{\Delta} have zero/non-zero structure (with non-zero coefficients given in the parameter \eqn{\Delta}). In case \code{diagCCtype = "beta13"} the diagonal elements are generated by a beta distribution with shape parameter equal to 1 and scale parameter equal to 3.
   
}
\value{
	An object of class \code{pcorSimJoint} containing the following components:

  \item{D1 }{dataset for first population.}
  \item{D2 }{dataset for second population.}
  \item{omega1 }{precision matrix for first population.}
  \item{omega2 }{precision matrix for second population.}
  \item{P }{total number of variables.}
  \item{diffs }{differential edges.}
  \item{delta }{generated values for the dependent structure.}
  \item{covJ }{joint covariance matrix used to generate the data.}
  \item{path1 }{adjacency matrix corresponding to the non-zero structure of \code{omega1}.}
  \item{path2 }{adjacency matrix corresponding to the non-zero structure of \code{omega2}.}
}
\references{
Cai, T., W. Liu, and X. Luo (2011). A Constrained L1 Minimization Approach to Sparse Precision Matrix Estimation. Journal of the American Statistical Association 106(494), 594-607. 

Newman, M. (2003). The structure and function of complex networks. SIAM REVIEW 45, 167-256. 

Wit, E. and A. Abbruzzo (2015, feb). Factorial graphical models for dynamic networks. Network Science 3(01), 37-57. 

Caballe, A., N. Bochkina, and C. Mayer (2016). Selection of the Regularization Parameter in Graphical Models using network charactaristics. eprint arXiv:1509.05326, 1-25. 
}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
} 
\seealso{
  \code{\link{pcorSimulator}}  for precision matrix generator. \cr
  \code{\link{plot.pcorSimJoint}}  for plotting joint partial correlation matrices.
}

\examples{
# example to use pcorSimulatorJoint function
EX1 <- pcorSimulatorJoint(nobs = 50, nclusters = 2, nnodesxcluster = c(30, 40), 
                          pattern = "pow", diffType = "cluster", dataDepend = "ind", 
                          pdiff = 0.2, diagCCtype = "dicot", diagNZ.strength = .5)
print(EX1)

EX2 <- pcorSimulatorJoint(nobs = 50, nclusters = 2, nnodesxcluster = c(30, 40), 
                          pattern = "pow", diffType = "rand", dataDepend = "diag", 
                          pdiff = 0.05, diagCCtype = "beta")
print(EX2)


}


