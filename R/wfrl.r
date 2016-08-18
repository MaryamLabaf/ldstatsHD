wfrl <- function(D1, D2, lambda1, lambda2, automLambdas=TRUE, paired = TRUE, sigmaEstimate = "CRmad", 
           maxiter=30, tol=1e-05, nsubset = 10000, rho = 1, rho.increment = 1, notOnlyLambda2 = TRUE)
{                           
     
     sigmaEstimate <- controlwfrl(D1 = D1, D2 = D2, lambda1 = lambda1, lambda2 = lambda2, paired = paired, 
    			automLambdas = automLambdas, sigmaEstimate = sigmaEstimate, maxiter = maxiter, 
    			tol = tol, nsubset = nsubset,  rho = rho, rho.increment = rho.increment, notOnlyLambda2=notOnlyLambda2)
	 
	 ## Initialisation 
     S		   <- list()
     S[[1]]    <- cor(D2[[1]])
     S[[2]]    <- cor(D2[[2]])
     P         <- dim(S[[1]])[2]
     N         <- dim(D1[[1]])[1]
     qp        <- dim(D1[[1]])[2]
     K         <- 2
     theta     <- list()
     thetae0   <- list()
     C12	   <- list()
 	 lengthPath <- length(lambda1)
 	     
     for (k in 1:K) {
        theta[[k]]  <- solve(S[[k]] + diag(rep(rho,P)))
        D1[[k]]     <- scale(D1[[k]])
        D2[[k]]     <- scale(D2[[k]])
        C12[[k]]    <- cor(D2[[k]],D1[[k]])
     }
    
	## initial iteration 
	if(automLambdas)  ind12 <- which(C12[[k]]!=1,arr.ind=TRUE)
    else              ind12 <- NULL
    if(!paired) corsOmega <- 0.5

    for (k in 1:K){
         thetae0[[k]]  <- theta[[k]] %*% (C12[[k]]*(N-1)/N)
    }
        
    if(paired){
		c12 			<- apply(as.matrix(1:P),1,function(d) theta[[1]][d,]/(theta[[1]][d,d]))
		c22 			<- apply(as.matrix(1:P),1,function(d) theta[[2]][d,]/(theta[[2]][d,d]))
		cors2 	    <- -apply(as.matrix(1:P),1,function(j) cor(D2[[1]][,j] + D2[[1]][,-j] %*% 
								c12[-j,j], D2[[2]][,j] + D2[[2]][,-j]%*%c22[-j,j]))

		cors1 	    <- -apply(as.matrix(1:qp),1,function(j) cor(D1[[1]][,j] - D2[[1]]%*%thetae0[[1]][,j], 
															   D1[[2]][,j] - D2[[2]]%*%thetae0[[2]][,j]))	    													   
        corsOmega 	<- t(cors1%*%t(cors2)) 
        
        #var1 	    <- apply(as.matrix(apply(as.matrix(1:qp),1,function(j) var(D1[[1]][,j] - D2[[1]]%*%thetae0[[1]][,j]))),1,function(x) rep(x,P))
        #var2 	    <- apply(as.matrix(apply(as.matrix(1:qp),1,function(j) var(D1[[2]][,j] - D2[[2]]%*%thetae0[[2]][,j]))),1,function(x) rep(x,P)) 
        
        #if(notOnlyLambda2) corsOmega 	<- t(cors1%*%t(cors2)) * (1/(sqrt(1 - (thetae0[[1]]*1)^2) * sqrt(1 - (thetae0[[2]]*1)^2)))
        
    }

	## calling recursive algorithm
    if(length(lambda1) > 1){
  	 objaux <- list()
  	 for(lamb1i in 1:lengthPath){  
     	pal           <- max(10^5,P*(P-1)/2)	   
        if(automLambdas & notOnlyLambda2& pal*lambda1[lamb1i]*lambda2>5)   
          lambda2 <- findingLambda2(pal, lam1=lambda1[lamb1i], lam2=lambda2, tol = 0.001, maxiter = 30) #internal use only
	 	
		objaux[[lamb1i]] <-   wfrlaux(D1 = D1, D2 = D2, lambda1 = lambda1[lamb1i], lambda2 = lambda2, automLambdas = automLambdas, 
		            paired = paired, sigmaEstimate = sigmaEstimate, maxiter = maxiter,  tol = tol, nsubset = nsubset, 
		            rho = rho, rho.increment = rho.increment, notOnlyLambda2 = notOnlyLambda2, corsOmega = corsOmega, 
		            thetae0 = thetae0, theta = theta, C12 = C12, ind12 = ind12)
	 }
	 obj 	    <- list()
	 obj$path   <- lapply(objaux, function(x) x$path)
	 obj$regCoef  <- lapply(objaux, function(x) x$regCoef)
	 obj$diff_value <- lapply(objaux, function(x) x$diff_value)
	 obj$iters <- lapply(objaux, function(x) x$iters)
	 rm(objaux)
	}
	else{
		pal           <- max(10^5,P*(P-1)/2)	   
        if(automLambdas & notOnlyLambda2& pal*lambda1*lambda2>5)   
          lambda2 <- findingLambda2(pal, lam1=lambda1, lam2=lambda2, tol = 0.001, maxiter = 30) #internal use only
		obj <-   wfrlaux(D1 = D1, D2 = D2, lambda1 = lambda1, lambda2 = lambda2, automLambdas = automLambdas, 
		            paired = paired, sigmaEstimate = sigmaEstimate, maxiter = maxiter,  tol = tol, nsubset = nsubset, 
		            rho = rho, rho.increment = rho.increment, notOnlyLambda2 = notOnlyLambda2, corsOmega = corsOmega, 
		            thetae0 = thetae0, theta = theta, C12 = C12, ind12 = ind12)
	}
	
	obj$paired  		<- paired
	obj$sigmaEstimate  	<- sigmaEstimate
    obj$lambda1  		<- lambda1
	obj$lambda2  		<- lambda2

    class(obj)  		<- "wfrl"
	return(obj)

}



wfrlaux <- function(D1, D2, lambda1, lambda2, automLambdas=FALSE, paired = FALSE, sigmaEstimate = "CRmad", 
           maxiter = 30, tol=1e-05, nsubset = 10000,  rho = 1, rho.increment = 1, notOnlyLambda2 = TRUE, 
           corsOmega = NULL, thetae0 = NULL, theta = NULL, C12 = NULL, ind12 = NULL)
{
           
	 ## Initialisation 
     P         <- dim(theta[[1]])[2]
     N         <- dim(D1[[1]])[1]
     qp        <- dim(D1[[1]])[2]
     K         <- 2
     thetae    <- thetae0
     Z         <- list()
     U         <- list()
     A         <- list()
     diff_value     <- numeric()
     diff_value[1]  <- 10
     thetae.prev1   <- matrix(0, P, qp)
     thetae.prev2   <- matrix(0, P, qp)
     
     iter 	   <- 0
     for (k in 1:K) {
        Z[[k]]      <- matrix(0, P, qp)
        U[[k]]      <- matrix(0, P, qp)
        A[[k]]      <- thetae[[k]] 
     }    
           
     ## finish initial iteration 
     Z       <- flsaFast(A, L = 1, lambda1, lambda2, corsOmega = corsOmega, P = P, ind12 = ind12, 
                       nsubset = nsubset, automLambdas = automLambdas, penalize.diagonal = TRUE,
                       sigmaEstimate = sigmaEstimate, FALSE, notOnlyLambda2 = notOnlyLambda2)
                      
	for (k in 1:K)    U[[k]]  <- U[[k]] + (thetae[[k]] - Z[[k]])
	
	iter               <-  iter + 1
	diff_value[iter]   <-  sum(abs(thetae[[1]] - thetae.prev1))/sum(abs(thetae.prev1)) +
						   sum(abs(thetae[[2]] - thetae.prev2))/sum(abs(thetae.prev2)) 
	thetae.prev1       <-  thetae[[1]]
	thetae.prev2       <-  thetae[[2]]
	rho                <-  rho * rho.increment
	
    ## recursive algorithm
    pb        	   <- txtProgressBar(min = 0, max = maxiter, style = 3)
    
    while (iter <= maxiter && diff_value[iter] > tol) {
        setTxtProgressBar(pb, iter)    
        
        for (k in 1:K){
           thetae[[k]]  <- theta[[k]] %*% (C12[[k]]*(N-1)/N  + Z[[k]] - U[[k]])  
           A[[k]]      <- thetae[[k]] + U[[k]]
        }

        Z       <- flsaFast(A, L = 1, lambda1, lambda2, corsOmega = corsOmega, P = P, ind12 = ind12, 
                       nsubset = nsubset, automLambdas = automLambdas, penalize.diagonal = TRUE,
                       sigmaEstimate = sigmaEstimate, FALSE, notOnlyLambda2 = notOnlyLambda2)
                      
        for (k in 1:K)    U[[k]]  <- U[[k]] + (thetae[[k]] - Z[[k]])
        
        iter               <-  iter + 1
        diff_value[iter]   <-  sum(abs(thetae[[1]] - thetae.prev1))/sum(abs(thetae.prev1)) +
                               sum(abs(thetae[[2]] - thetae.prev2))/sum(abs(thetae.prev2)) 
        thetae.prev1       <-  thetae[[1]]
        thetae.prev2       <-  thetae[[2]]
        rho                <-  rho * rho.increment
    }
    close(pb)
	
	path 		<- list()
	path[[1]] 	<- (Z[[1]]!=0)*1
	path[[2]] 	<- (Z[[2]]!=0)*1

    obj 				<- list(regCoef = Z, path = path, diff_value = diff_value, iters = iter, corsOmega = corsOmega)
	
    return(obj)
}


