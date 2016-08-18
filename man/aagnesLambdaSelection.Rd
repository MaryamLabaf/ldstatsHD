\name{amseLambdaSelection}
\alias{amseLambdaSelection}


\title{
  Augmented-MSE regularization parameter selection
 }
\description{

	\code{amseLambdaSelection} is a function designed to select the regularization 
  parameter lambda in graphical models that compromises global clustering structure 
  and variability of the graph.
}
\usage{

amseLambdaSelection(obj, pathIni, y, generator = c("subsampling", "montecarlo"),
                    pB = 0.7, nite = 10, method = "mb", from = 1, until = NULL, 
                    distF = c("correlation","shortPath"), oneByone = FALSE, 
                    many = 3)
}

\arguments{

\item{obj}{
      an object of class \code{\link{huge}} or \code{\link{camel.tiger}}.
}
\item{pathIni}{
      path with best global characteristics (for instance the \code{\link{agnesLambdaSelection}} selected path).
}
\item{y}{
      original \eqn{n\times p} data set.
}
\item{generator}{
	type of generator to find the mean squared error: name that uniquely identifies  \code{"subsampling"} or \code{"montecarlo"}.
}
\item{pB}{
      proportion of observations used in subsampling iterations.
}
\item{nite}{
      number of iterations used to estimate the mean square error.
}
\item{method}{
      method used to estimate the networks: name that uniquely identifies \code{"mb"}, \code{"glasso"} or \code{"tiger"}. 
}

\item{from}{
      starting point in lambda sequence. 
}
\item{until}{
      last point in lambda sequence. If \code{until = NULL}, all lambda sequence is explored.
}

\item{distF}{
      distance function used to find the dissimilarity matrix from the graph: name that uniquely identifies \code{"correlation"} or \code{"shortPath"}.
}
\item{oneByone}{
      If \code{TRUE}, the estimation process is done separately for each \eqn{\lambda}.
}

\item{many}{
      If \code{oneByone = TRUE}, the estimation process is done separately for every \code{many} \eqn{\lambda}'s. 
}


}
\details{
A-MSE algorithm finds \eqn{\lambda} by minimizing the risk function
\deqn{
R_{AMSE}(\lambda) = E(\sum_{i>j} |\delta_{ij}-\hat{\delta}_{ij}^{\lambda}|^2) 
}
where \eqn{\hat{\delta}_{ij}^{\lambda}} is the dissimilarity matrix of the graph 
(see \code{\link{graphCorr}}). The expected value is approximated by either subsampling or 
Monte Carlo and the theoretical \eqn{\delta_{ij}} is approximated by the graph in \code{pathIni}. 

We recommend using the AGNES algorithm \code{\link{agnesLambdaSelection}} to 
approximate \code{pathIni} since provides good reference of global network 
structure for clustered-based graph structures. Then, A-MSE gives a good trade-off 
between graph variability and global network characteristics.

If \code{pathIni} is given  by \code{\link{agnesLambdaSelection}} and 
\code{generator = "subsampling"}, then the lambda selected is always smaller than 
the lambda obtained by \code{AGNES}.
      
The \code{oneByone} approach is suggested to save memory space for very high-dimensional data. 
}
\value{          
	An object of class \code{lambdaSelection} containing the following components: 
 
  \item{opt.lambda }{optimal lambda.} 
  \item{crit.coef }{coefficients for each lambda given the criterion A-MSE.}
  \item{criterion }{with value \code{"A-MSE"}.}
  

}
\references{
Caballe, A., N. Bochkina, and C. Mayer (2016). Selection of the Regularization Parameter in Graphical Models using network charactaristics. eprint arXiv:1509.05326, 1-25. 
}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}
\seealso{
 \code{\link{lambdaSelection}} for other lambda selection approaches.\cr
}

\examples{

# example to use amse function
EX1         <- pcorSimulator(nobs = 70, nclusters = 3, nnodesxcluster = c(40,30,20), 
                             pattern = "powerLaw")
y           <- EX1$y
Lambda.SEQ  <- seq(.25, 0.70, length.out = 40)
out3        <- huge(y, method = "mb", lambda = Lambda.SEQ)
AG.COEF     <- agnesLambdaSelection(out3, distF = "shortPath", way = "direct")
AG.LAMB     <- which(AG.COEF$opt.lambda == Lambda.SEQ)

## not run
#AAG.COEF    <- amseLambdaSelection(out3, out3$path[[AG.LAMB]], y = y, 
#                 distF = "shortPath", from = AG.LAMB)
#print(AAG.COEF) 

}


