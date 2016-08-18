\name{cor2mean.adj}
\alias{cor2mean.adj}


\title{
  adjusted average square correlation by rows
 }
\description{
 Finds in a computationally fast algorithm the adjusted average square correlation magnitude for every variable of a dataset.  
}
\usage{

cor2mean.adj(mat)

}

\arguments{

\item{mat}{
      \eqn{p \times n} matrix with the p-variate dataset.
}

}
\details{
 The adjusted average square correlation of variable \eqn{i} is given by 
 \deqn{
  (n-1)/(n-2) \bar{r}_{i}^2 - 1/(n-2)
 }
 where \eqn{n} is the sample size and \eqn{\bar{r}_{i}^2} is the average square 
 correlation matrix for the \eqn{i}th row, which is computed by \code{\link{cor2mean}}.
 
}
\value{

 A vector containing the adjusted square average correlation (excluding diagonal) for every variable.

}

\references{
To come
}
\author{
 Mayer, Claus, Adria Caballe and Natalia Bochkina.
}
\seealso{
 \code{\link{cor2mean}} for average square correlations.}

\examples{
EX1        <- pcorSimulator(nobs = 50, nclusters= 3, nnodesxcluster = c(100,30,50), 
                            pattern = "powerLaw", plus = 0)
corsEX1     <- cor2mean(t(EX1$y))
corsadjEX1  <- cor2mean.adj(t(EX1$y))

}
