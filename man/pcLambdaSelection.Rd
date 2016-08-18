\name{pcLambdaSelection}
\alias{pcLambdaSelection}

\title{
  Path Connectivity regularization parameter selection
 }
\description{

	\code{pcLambdaSelection} is a function designed to select the regularization 
  parameter in graphical models. It selects the graph which captures the biggest drop
  in graph connectivity.
	
}
\usage{

pcLambdaSelection(obj)

}

\arguments{

\item{obj}{
      an object of class \code{\link{huge}}, \code{\link{camel.tiger}} or \code{\link{wfgl}}. 
}
}
\details{

Path Connectivity (PC) algorithm finds \eqn{\lambda}  by maximizing the biggest drop of 
connectivity in estimated graphs. We define connectivity by the average geodesic distance between pairs of nodes
(see \code{\link{graphDist}}). 

PC gives a fast and suitable way to select \eqn{\lambda}  when there are distinct 
clusters in the data. Given two graphs, corresponding to two consecutive \eqn{\lambda}'s, 
the difference between the average geodesic distance will be large if the first graph 
contains edges that connect different clusters which are not present in the second graph.

Note that PC should  be used when fitting graphical models with an equidistant sequence
for \eqn{\lambda}.



}
\value{          
	An object of class \code{lambdaSelection} containing the following components: 
 
  \item{opt.lambda }{optimal lambda.} 
  \item{crit.coef }{coefficients for each lambda given the criterion PC.}
  \item{criterion }{with value \code{"PC"}.}
  

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

# example to use pc function
EX1         <- pcorSimulator(nobs = 70, nclusters = 3, nnodesxcluster = c(40,30,20), 
                             pattern = "powerLaw")
y           <- EX1$y
Lambda.SEQ  <- seq(.25,0.70,length.out = 40)
out3        <- huge(y, method = "mb", lambda = Lambda.SEQ)
PC.COEF     <- pcLambdaSelection(out3)
print(PC.COEF)
 
 
}


