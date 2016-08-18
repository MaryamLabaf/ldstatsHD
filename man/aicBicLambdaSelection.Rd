\name{aicAndbicLambdaSelection}
\alias{aicAndbicLambdaSelection}


\title{
  AIC/BIC regularization parameter selection
 }
\description{

	\code{aicAndbicLambdaSelection} is a function designed to select the regularization 
  parameter in graphical models. It selects the graph with smallest AIC or BIC coefficients.
	
}
\usage{

aicAndbicLambdaSelection(obj,  y, criterion = c("AIC", "BIC"))
}

\arguments{

\item{obj}{
      an object of class \code{\link{huge}} or \code{\link{camel.tiger}}. 
}
\item{y}{
      original \eqn{n \times p} data set.
}
\item{criterion}{
	coefficients and optimal lambdas to be stored: to select from \code{"AIC"}, \code{"BIC"} or both.
}

}
\value{          
	An object of class \code{lambdaSelection} containing the following components: 
 
  \item{opt.lambda }{optimal lambdas for AIC and BIC.} 
  \item{crit.coef }{coefficients for each lambda given the criterion AIC and BIC.}
  \item{criterion }{with value defined by argument \code{criterion}.}
  

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

# example to use aicAndBic function
EX1         <- pcorSimulator(nobs = 50, nclusters = 3, nnodesxcluster = c(40,30,30), 
                             pattern = "powerLaw")
y           <- EX1$y
Lambda.SEQ  <- seq(.35, 0.70, length.out = 40)
out3        <- huge(y, method = "glasso", lambda = Lambda.SEQ, cov.output = TRUE)
AIC.COEF    <- aicAndbicLambdaSelection(out3, y = y)
print(AIC.COEF)

 
 
}

