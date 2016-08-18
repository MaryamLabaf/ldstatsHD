\name{ldstatsHD-package}
\alias{ldstatsHD-package}
\alias{ldstatsHD}
\docType{package}

\title{
  Linear Dependence Statistics for High-Dimensional data
}
\description{
    This package consists of functions with statistical methods related 
	to the estimation and testing of multiple correlation, partial correlation and
	regression coefficient matrices when data is high-dimensional (\eqn{n<p}).	
   
	Joint estimation of two partial correlation matrices (see \code{\link[ldstatsHD]{wfgl}}) 
	and joint estimation of two regression coefficient matrices (see \code{\link[ldstatsHD]{wfrl}}) 
	are currently implemented in the package. These use a weighted-fused lasso penalized maximum 
	likelihood estimator such that they encourage both sparsity and similarity between estimated
	matrices. 
	
	\pkg{ldstatsHD} also contains approaches to select the sparsity tuning parameter of graphical 
	lasso estimators such that several risk functions based on characteristics of the estimated 
	networks are available (see \code{\link[ldstatsHD]{lambdaSelection}}). Among others, statistics that
	measure clustering structure or network connectivity can be used to find the desired networks.
	
	It finally includes statistical methods that test global dependence characteristics:
	(i) a test for equality of two correlation matrices as well as a test for Identity correlation matrix
	(see \code{\link[ldstatsHD]{eqCorrMatTest}}); (ii)  a test for equality of two correlation matrix rows
	as well as a test to check if a variable is linearly independent of the rest of the variables in a dataset
	(see \code{\link[ldstatsHD]{eqCorTestByRows}}).
	
	A particularity of the implemented methods in \pkg{ldstatsHD} is that it permits cases where datasets are dependent (e.g. paired data). 

	\pkg{ldstatsHD} provides two partial correlation matrix simulators such that all methods can be tested using  
	 using simulated data: see  \code{\link[ldstatsHD]{pcorSimulator}} to generate a single partial correlation / dataset
	 and \code{\link[ldstatsHD]{pcorSimulatorJoint}} to generate a joint partial correlation matrix and two (dependent) datasets.
}
\details{
\tabular{ll}{
Package: \tab ldstatsHD\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2016-07-08\cr
License: \tab 	GPL-2\cr
LazyLoad: \tab yes\cr     
}
}
\author{
 Caballe, Adria <a.caballe@sms.ed.ac.uk>, Natalia Bochkina and Claus Mayer.
}
\references{
To come

}



