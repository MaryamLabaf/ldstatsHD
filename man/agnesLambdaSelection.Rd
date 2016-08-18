\name{agnesLambdaSelection}
\alias{agnesLambdaSelection}


\title{
  AGNES regularization parameter selection
 }
\description{

	\code{agnesLambdaSelection} is a function designed to select the regularization 
  parameter in graphical models. It selects the most clustered conditional dependence graph structure where clusters are
  defined by the hierarchical algorithm \code{\link{agnes}} (See details).
	
}
\usage{

agnesLambdaSelection(obj, way = "direct", nite = 10, subsvec = NULL,
                     eps = 0.05, until = NULL, minNodes = 30, 
                     distF = c("correlation","shortPath"))
                     
}
\arguments{

\item{obj}{
      an object of class \code{\link{huge}}, \code{\link{camel.tiger}} or \code{\link{wfgl}}. 
}
\item{way}{
      name that uniquely identifies \code{"direct"} (default), \code{"rand.sampling"} for random subsets algorithm 
      and \code{"int.sampling"} for intelligent subsets algorithm.
}
\item{nite}{
      vector with the number of iterations used for each lambda  (only if \code{way = "rand.sampling"} or \code{way = "int.sampling"}).
}
\item{subsvec}{
      vector with the number of subsamples used for each lambda (only if \code{way = "rand.sampling"} or \code{way = "int.sampling"}). If \code{NULL}, 
      argument \code{minNodes} determines the number of subsamples for all lambdas. 
}
\item{eps}{
      acceptance tolerance for subsets of variables.
}
\item{until}{
      the last path used in \code{obj}. If \code{NULL}, all paths are used to select lambda. 
}
\item{minNodes}{
      minimum number of nodes with connections to compute the AGNES coefficient 
      (the coefficient is zero for paths with less nodes than \code{minNodes}).
}

\item{distF}{
      distance function used to find the dissimilarity matrix from the graph: name that uniquely identifies 
      \code{"correlation"} and \code{"shortPath"}.
}

}
\details{
AGNES algorithm finds \eqn{\lambda}  by minimizing the risk function
\deqn{
R_{AGNES}(\lambda) =  -AC(\lambda)
}
where  \eqn{AC(\lambda)} is the AGNES coefficient calculated using the R function 
\code{\link{agnes}}. Using AGNES we select the \eqn{\lambda} that maximizes the between 
vs within cluster dissimilarities ratio given the dissimilarity matrix of the graph 
(see \code{\link{graphCorr}} and \code{\link{graphDist}} for possible dissimilarities).

A variable subset selection algorithm is available to estimate \eqn{AC(\lambda)} for 
very high-dimensional data. It is recommended in order to save memory space and computational time.
Especially \code{way = "int.sampling"} which tends to finds similar lambda selections to the default 
procedure.

\code{agnesLambdaSelection} gives a good recovery of global network characteristics when the true partial correlation matrix is block diagonal.
}
\value{          
	An object of class \code{lambdaSelection} containing the following components: 
 
  \item{opt.lambda }{optimal lambda.} 
  \item{crit.coef }{coefficients for each lambda given the criterion AGNES.}
  \item{criterion }{with value \code{"AGNES"}.}
  

}
\references{
Caballe, A., N. Bochkina, and C. Mayer (2016). Selection of the Regularization Parameter in Graphical Models using network charactaristics. eprint arXiv:1509.05326, 1-25. 
}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}
\seealso{
 \code{\link{lambdaSelection}} for other lambda selection approaches and \code{\link{agnes}}
 for clustering implementation.
}

\examples{

# example to use agnes function
EX1         <- pcorSimulator(nobs = 70, nclusters = 3, nnodesxcluster = c(40,30,20), 
                             pattern = "powerLaw")
y           <- EX1$y
Lambda.SEQ  <- seq(.25, 0.70, length.out=40)
out3        <- huge(y, method = "mb", lambda = Lambda.SEQ)
AG.COEF     <- agnesLambdaSelection(out3, distF = "shortPath", way = "direct")
print(AG.COEF)
 
 
}

