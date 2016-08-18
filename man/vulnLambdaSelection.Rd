\name{vulLambdaSelection}
\alias{vulLambdaSelection}


\title{
  Vulnerability regularization parameter selection
 }
\description{

	\code{vulLambdaSelection} is a function designed to select the regularization 
  parameter in graphical models. It selects the graph with largest average nodes vulnerability.
	
}
\usage{

vulLambdaSelection(obj, loo = FALSE, subOut =  10, nite = 50)

}

\arguments{

\item{obj}{
      an object of class \code{\link{huge}}, \code{\link{camel.tiger}} or \code{\link{wfgl}}. 
}
\item{loo}{
      if \code{TRUE} an exhaustive leave-one-out procedure is done, 
      otherwise it is used a subsampling approach with \code{nite} iterations and 
      leaving out \code{subOut} variables. 
}
\item{subOut}{
      number of variables left out in each iteration (only used if \code{loo = FALSE}).
}
\item{nite}{
      number of iterations (only used if \code{loo = FALSE}).
}


}
\details{
Vulnerability algorithm finds lambda by minimizing the risk function
\deqn{
R_{VUL}(\lambda) =  - \sum_{i=1}^p \frac{E^\lambda-E_i^\lambda}{E^\lambda}
}
where  \eqn{E^\lambda} is the global efficiency of the original network and
\eqn{E^\lambda_i} is the global efficiency of the network once eliminating the node
\eqn{i}. Global efficiency is defined by the harmonic mean of the geodesic distance 
 (see \code{\link{graphDist}}). 

Vulnerability gives \eqn{\lambda} selection that contains the most vulnerable graph,
meaning that the removal of a node in the network in average would affect the most the 
estimated graph.
}
\value{          
  An object of class \code{lambdaSelection} containing the following components: 
 
  \item{opt.lambda }{optimal lambda.} 
  \item{crit.coef }{coefficients for each lambda given the criterion VUL.}
  \item{criterion }{with value \code{"VUL"}.}
  
}
\references{
Costa, L. and F. Rodrigues (2007). Characterization of complex networks: A survey of measurements. Advances in Physics 56(1), 167-242. 
}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}
\seealso{
 \code{\link{lambdaSelection}} for other lambda selection approaches.\cr
}

\examples{

# example to use vul function
EX1         <- pcorSimulator(nobs = 50, nclusters = 2, nnodesxcluster = c(40,30), 
                             pattern="powerLaw")
y           <- EX1$y
Lambda.SEQ  <- seq(.35, 0.70, length.out = 10)
out3        <- huge(y, method = "mb", lambda = Lambda.SEQ)
## not run
#VUL.COEF    <- vulLambdaSelection(out3)
#print(VUL.COEF)

 
 
}


