controlwfrl <- function(D1, D2, lambda1, lambda2, paired, automLambdas, 
	 			 sigmaEstimate, maxiter, tol, nsubset,
	 			 rho, rho.increment,notOnlyLambda2)
{
  ## dimensions
  if(dim(D1[[1]])[2] != dim(D1[[2]])[2])
    stop("Response datasets dimensions must be equal")

  if(dim(D2[[1]])[2] != dim(D2[[2]])[2])
    stop("Covariates datasets dimensions must be equal")
  
  if(dim(D1[[1]])[1] != dim(D2[[1]])[1])
    stop("Response and covariates must have same number of observations (pop.1)")

  if(dim(D1[[2]])[1] != dim(D2[[2]])[1])
    stop("Response and covariates must have same number of observations (pop.2)")
  
  if(dim(D1[[1]])[1] != dim(D1[[2]])[1] & paired)
    stop("Response sample sizes must be equal when paired = TRUE")

  if(dim(D1[[2]])[1] != dim(D2[[2]])[1] & paired)
    stop("covariates sample sizes must be equal when paired = TRUE")
  
  
  ## lambda1, lambda2, automLambdas
  if(length(lambda2) > 1) 
    stop("lambda2 is not well defined. Must be a single value")
  if( any(lambda1 < 0) & !automLambdas)
    stop("lambda1 is not well defined for automLambda = FALSE. It must be larger or equal to zero")
  if( lambda2 < 0 & !automLambdas)
    stop("lambda2 is not well defined for automLambda = FALSE. It must be larger or equal to zero")  
  if( any(lambda1 < 0 | lambda1 > 0.5 ) & automLambdas &notOnlyLambda2)
    stop("lambda1 is not well defined for automLambda = TRUE. It must be between 0 and 0.5")
  if( (lambda2 < 0 | lambda2 > 0.5 ) & automLambdas)
    stop("lambda2 is not well defined for automLambda = TRUE. It must be between 0 and 0.5")    	
  
  ## prospective sigmaEstimate
  if(automLambdas){
   psigmaEstimate <- c("CRmad", "mad", "IQR")
   if(length(sigmaEstimate) > 1)
    warning("sigmaEstimate attribute length is larger than one. Only the first component will be used")
 
   sigmaEstimate <- psigmaEstimate[pmatch(sigmaEstimate[1], psigmaEstimate)]
   if (is.na(sigmaEstimate)) stop("sigmaEstimate is not well define. It must be selected from \"CRmad\", \"mad\" or \"IQR\" ")
  }
  else
   sigmaEstimate <- "CRmad"
   
 
  ## maxiter, tol, nsubset, rho, rho.increment
  if(maxiter < 1)
    stop("not enought permutations. Increase nite to be larger than 1")
    
  if(tol <= 0)
    stop("tol must be larger than zero")

  if(!is.null(nsubset) & automLambdas)
  {
   if(nsubset <= 50)
     stop("nsubset is very small. It should be larger than 50")
  }
  
  if(rho.increment < 0)
     stop("rho.increment must be larger than zero")	

  return(c(sigmaEstimate))
}
