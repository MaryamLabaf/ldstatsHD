\name{pcorSimulator}
\alias{pcorSimulator}
\alias{print.pcorSim}
\alias{summary.pcorSim}


\title{
  Partial Correlation Matrix simulator
 }
\description{

 \code{pcorSimulator} creates a block diagonal positive definite precision matrix with three possible 
   graph structures: hubs-based, power-law and random. Then, it generates samples from a multivariate normal 
   distribution with covariance matrix given by the inverse of such precision matrix.
   
}
\usage{

pcorSimulator(nobs, nclusters, nnodesxcluster, pattern = "powerLaw", 
              low.strength = 0.5, sup.strength = 0.9, nhubs = 5, 
              degree.hubs = 20, nOtherEdges = 30, alpha = 2.3, plus = 0, 
              prob = 0.05, perturb.clust = 0, mu = 0,
              probSign = 0.5, seed = sample(10000, nclusters))
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
      graph structure pattern: name that uniquely identifies \code{"hubs"}, \code{"powerLaw"} and \code{"random"}.
}
\item{low.strength}{
      minimum magnitude for nonzero partial correlation elements before regularization.
}
\item{sup.strength}{
      maximum magnitude for nonzero partial correlation elements before regularization.
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
      probability of edge presence for random networks (if \code{pattern = "random"}).
}                                                                                            
\item{perturb.clust}{
     proportion of the total number of edges that are connecting two different clusters.
}
\item{mu}{
     expected values vector to generate data (zero by default).
}
\item{probSign}{
     probability of positive sign for non-zero partial correlation coefficients. Thus, negative signs 
     are obtained with probability \code{1-probSign}. 
}
\item{seed}{
      vector with seeds for each cluster.
}
}
\details{
Hubs-based networks are graphs where only few nodes have a much higher degree 
(or connectivity) than the rest. Power-law networks assume that the variable \eqn{p_k}, 
which denotes the fraction of nodes in the network that has degree \eqn{k},
is given by a power-law distribution
\deqn{
p_k =  \frac{k^{-\alpha}}{\varsigma(\alpha)},
}
for \eqn{k \geq 1}, a constant \eqn{\alpha>0} and the normalizing function \eqn{\varsigma(\alpha)} 
which is the Riemann zeta function. Finally, random networks are also defined by the distribution in 
the proportion \eqn{p_k}. In this case, \eqn{p_k} follows a binomial distribution
\deqn{
p_k =  {p\choose k} \theta^k (1-\theta)^{p-k},
}
where the parameter \eqn{\theta} determines the proportion of edges (or sparsity) in the graph. 

The regularization is given by \eqn{\Omega^{(1)} = \Omega^{(0)} +  \delta I}, with \eqn{\delta}
such that the condition number of  \eqn{\Omega^{(1)}} is less than the number of nodes.

}
\value{       

	An object of class \code{pcorSim} containing the following components: 
 
  \item{y }{generated data set.} 
  \item{hubs }{hub nodes position.}
  \item{edgesInGraph }{edges given by the non-zero elements in the precision matrix.}
  \item{omega }{precision matrix used to generate the data.}
  \item{covMat }{covariance matrix used to generate the data.}
  \item{path }{adjacency matrix corresponding to the non-zero structure of \code{omega}.}
  
              
}
\references{
Cai, T., W. Liu, and X. Luo (2011). A Constrained L1 Minimization Approach to Sparse Precision Matrix Estimation. Journal of the American Statistical Association 106(494), 594-607. 

Newman, M. (2003). The structure and function of complex networks. SIAM REVIEW 45, 167-256. 

Caballe, A., N. Bochkina, and C. Mayer (2016). Selection of the Regularization Parameter in Graphical Models using network charactaristics. eprint arXiv:1509.05326, 1-25. 
}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}

\seealso{
\code{\link{plot.pcorSim}} for graphical representation of the generated partial correlation matrix.\cr
\code{\link{pcorSimulatorJoint}} for joint partial correlation matrix generation.
}
\examples{
# example to use pcorSimulator function

EX1 <- pcorSimulator(nobs = 50, nclusters=3, nnodesxcluster=c(100,30,50), 
                pattern="powerLaw", plus=0)
print(EX1)
                
EX2 <- pcorSimulator(nobs = 25, nclusters=2, nnodesxcluster=c(60,40), 
                pattern="powerLaw", plus=1)
print(EX2)
 
}


