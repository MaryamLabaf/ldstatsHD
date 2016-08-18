\name{lambdaSelection}
\alias{lambdaSelection}
\alias{print.lambdaSelection}
\alias{summary.lambdaSelection}
\alias{plot.lambdaSelection}


\title{
  Regularization parameter selection based on network characteristics
 }
\description{

	\code{lambdaSelection} is a function designed to select the sparsity regularization 
  parameter in graphical models.
	
  There are available seven different criterion to select lambda with risk functions
  based on network characteristics: Path connectivity (PC), AGlommerative NESted (AGNES), Augmented-MSE (A-MSE),
  Vulnerability (VUL), AIC/BIC and StARS (from \code{\link{huge}} package).
 
}
\usage{

lambdaSelection(obj, criterion = c("PC","AGNES","A-MSE","VUL","STARS", "AIC", "BIC"), 
                ...)

}

\arguments{

\item{obj}{
      an object of class \code{\link{huge}}, \code{\link{camel.tiger}} or \code{\link{wfgl}}. 
}
 \item{criterion}{
  regularization parameter selection approach: name that uniquely identifies \code{"PC"} (Path connectivity), 
\code{"AGNES"} (AGlommerative NESted Algorithm), \code{"A-MSE"} (Augmented mean square error), 
\code{"VUL"} (Maximum VULnerability), \code{"AIC"}/\code{"BIC"} (minimum AIC/BIC) or \code{"StARS"} (Stability approach).
}

  \item{\dots}{
    arguments passed to or from other methods to the low level. See \code{\link{pcLambdaSelection}},
\code{\link{agnesLambdaSelection}}, \code{\link{amseLambdaSelection}}, \code{\link{vulLambdaSelection}},\cr
\code{\link{aicAndbicLambdaSelection}} and \code{\link{huge}} for details.
  }   
}
\details{
This function considers seven ways of selecting the regularization parameter in 
graphical models by minimizing a certain risk function based only on network characteristics
of the underlying structure of \eqn{\Omega}

\deqn{
\hat\lambda = \arg \min_{\lambda} R(\lambda, \hat{G}_\lambda),
}
where \eqn{\hat{G}_\lambda} is the estimated graph structure of \eqn{\hat{\Omega}}. For 
instance see \code{\link{pcLambdaSelection}}, \cr
\code{\link{agnesLambdaSelection}}, \code{\link{amseLambdaSelection}}, \code{\link{vulLambdaSelection}}, \cr
\code{\link{aicAndbicLambdaSelection}} and \code{\link{huge}} for the implemented  criterions to select \eqn{\lambda}.

For \code{\link{wfgl}} objects, only \code{criterion = c("PC","AGNES","VUL")} are available.

}
\value{          
  An object of class \code{lambdaSelection} containing the following components: 
  
  \item{opt.lambda }{optimal lambda.} 
  \item{crit.coef }{coefficients for each lambda given the criterion.}
  \item{criterion }{criterion used to select lambda.}
  
}
\references{
Caballe, A., N. Bochkina, and C. Mayer (2016). Selection of the Regularization Parameter in Graphical Models using network charactaristics. eprint arXiv:1509.05326, 1-25. 
}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}
\note{
For large dimensions, \code{"A-MSE"}, \code{"VUL"} (the most) and \code{"StARS"} can be computationally intensive. 
}
\seealso{
 \code{\link{pcLambdaSelection}}, \code{\link{agnesLambdaSelection}}, \code{\link{amseLambdaSelection}}\cr
 \code{\link{vulLambdaSelection}}, \code{\link{aicAndbicLambdaSelection}} and \code{\link{huge}}.\cr
}

\examples{

# example to use agnes function
EX1         <- pcorSimulator(nobs = 50, nclusters = 3, nnodesxcluster = c(40,30,30), 
                             pattern="powerLaw")
y           <- EX1$y
Lambda.SEQ  <- seq(.35,0.70, length.out = 40)
out3        <- huge(y, method = "mb", lambda = Lambda.SEQ)
PC.COEF     <- lambdaSelection(out3, criterion = "PC")
#AG.COEF     <- lambdaSelection(out3, criterion = "AGNES")
 
}


