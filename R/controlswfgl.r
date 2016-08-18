controlwfgl <- function(D1, D2, lambda1, lambda2, paired, automLambdas, 
	 			 sigmaEstimate, pairedEst, maxiter, tol, nsubset, weights,
	 			 rho, rho.increment, triangleCorrection, alphaTri, temporalFolders)
{
  ## dimensions
  if(dim(D1)[2] != dim(D2)[2])
    stop("Datasets dimensions must be equal")
  
  if(dim(D1)[1] != dim(D2)[1] & paired)
    stop("Datasets sample sizes must be equal when paired = TRUE")
  
  ## lambda1, lambda2, automLambdas
  if(length(lambda2) > 1) 
    stop("lambda2 is not well defined. Must be a single value")
  if( any(lambda1 < 0) & !automLambdas)
    stop("lambda1 is not well defined for automLambda = FALSE. It must be larger or equal to zero")
  if( lambda2 < 0 & !automLambdas)
    stop("lambda2 is not well defined for automLambda = FALSE. It must be larger or equal to zero")  
  if( any(lambda1 < 0 | lambda1 > 0.5 ) & automLambdas)
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
   
  ## prospective pairedEst
  if(paired){
   ppairedEst <- c("Reg-based-sim", "Reg-based")
   if(length(pairedEst) > 1)
    warning("pairedEst attribute length is larger than one. Only the first component will be used")
 
   pairedEst <- ppairedEst[pmatch(pairedEst[1], ppairedEst)]
   if (is.na(pairedEst)) stop("pairedEst is not well define. It must be selected from \"Reg-based-sim\" or \"Reg-based\" ")
  }
  else
   pairedEst <- "Reg-based-sim"
 
  ## maxiter, tol, nsubset, weights, rho, rho.increment
  if(maxiter < 1)
    stop("not enought iterations. Increase maxiter to be larger than 1")
    
  if(tol <= 0)
    stop("tol must be larger than zero")

  if(!is.null(nsubset) & automLambdas)
  {
   if(nsubset <= 50)
     stop("nsubset is very small. It should be larger than 50")
  }
  
  if(length(weights)!=2)
     stop("weights must be a vector of length 2")
     
  if(any(weights <= 0))
     stop("weights must be larger than zero")
     
  if(rho.increment < 0)
     stop("rho.increment must be larger than zero")
   
  ## triangleCorrection, alphaTri, temporalFolders
  if(triangleCorrection){
   if( alphaTri < 0 | alphaTri > 0.5 ) 
     stop("alphaTri is not well defined for triangleCorrection = TRUE. It must be between 0 and 0.5")
  }  	

  return(c(sigmaEstimate, pairedEst))
}
