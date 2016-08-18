\name{cor2mean}
\alias{cor2mean}

\title{
  Average square correlation by rows
 }
 
\description{
  Finds in a computationally fast algorithm the average square correlation magnitude for every variable of a dataset.  
  }
\usage{
cor2mean(mat)
}

\arguments{

\item{mat}{
      \eqn{p \times n} matrix with the p-variate dataset.
}
}
\details{
	  It is especially suitable for high dimensions. For instance it handles well dimensions of order of thousands.
}

\value{          
 The average square correlation magnitude of the sample correlation matrix (including the diagonal) for every variable in \code{mat}.
}
\references{
To come
}
\author{
 Mayer, Claus, Adria Caballe and Natalia Bochkina.
}

\seealso{
  \code{\link{cor2mean.adj}} for adjusted average square correlation magnitude.
}
\examples{
EX1      <- pcorSimulator(nobs = 50, nclusters= 3, nnodesxcluster = c(100,30,50), 
                          pattern = "powerLaw", plus = 0)
corsEX1  <- cor2mean(t(EX1$y))

}


