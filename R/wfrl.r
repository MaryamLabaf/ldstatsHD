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
    if(!paired) corsOmega <- array(0,dim=c(P,qp))

    for (k in 1:K){
         thetae0[[k]]  <- theta[[k]] %*% (C12[[k]]*(N-1)/N)
    }
        
    if(paired){
		c12 			<- apply(as.matrix(1:P),1,function(d) theta[[1]][d,]/(theta[[1]][d,d]))
		c22 			<- apply(as.matrix(1:P),1,function(d) theta[[2]][d,]/(theta[[2]][d,d]))
		cors2 	        <- -apply(as.matrix(1:P),1,function(j) cor(D2[[1]][,j] + D2[[1]][,-j] %*% 
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
		obj <-   wfrlaux(D1 = D1, D2 = D2, lambda1 = lambda1, lambda2 = lambda2, automLambdas = automLambdas, 
		            paired = paired, sigmaEstimate = sigmaEstimate, maxiter = maxiter,  tol = tol, nsubset = nsubset, 
		            rho = rho, rho.increment = rho.increment, notOnlyLambda2 = notOnlyLambda2, corsOmega = corsOmega, 
		            thetae0 = thetae0, theta = theta, C12 = C12, ind12 = ind12)
	}
	
	obj$paired  		<- paired
	obj$sigmaEstimate  	<- sigmaEstimate
    obj$lambda1  		<- lambda1
	obj$lambda2  		<- lambda2
    obj$automLambdas	<- automLambdas	

    class(obj)  		<- "wfrl"
	return(obj)

}



wfrlaux <- function(D1, D2, lambda1, lambda2, automLambdas=FALSE, paired = FALSE, sigmaEstimate = "CRmad", 
           maxiter = 30, tol=1e-05, nsubset = 10000,  rho = 1, rho.increment = 1, notOnlyLambda2 = TRUE, 
           corsOmega = NULL, thetae0 = NULL, theta = NULL, C12 = NULL, ind12 = NULL )
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
           
     
    if(lambda2>0&automLambdas){
     sq1  <- seq(-8,8,length.out=10000)
	 alpha2 <- lambda2
	 thetaseq <- seq(-0.3,0.99,length.out=1000)
		fas <- apply(as.matrix(thetaseq),1,function(theta){
         x2 <- sort(abs(sq1-qnorm(1-alpha2/2)*sqrt(2-2*theta)/2),index.return=TRUE)
		 x2$x[which(cumsum( (dnorm(sq1)*diff(sq1)[1] * pnorm((sq1*(1-theta) -qnorm(1-alpha2/2)*sqrt(2-2*theta))/sqrt(1-theta^2))/(alpha2/2))[x2$ix])>1-lambda1)[1]]
     })
		iddd 		<- pmax(1,round(punif(corsOmega,-0.3,0.99)*1000))
		alpha1A 	<- matrix(fas[iddd],nrow=dim(corsOmega)[1]) 
		
		thetaSeq2 <- sample(1:length(corsOmega),1000)
 		fas2 <- apply(as.matrix(thetaSeq2),1,function(j){
 			theta <- corsOmega[j]
 			lam1  <- alpha1A[j]
			sq1  <- seq(lam1 + qnorm(1-alpha2/2)*sqrt(2-2*theta)/2,10,length.out=10000)
			da1  <- dnorm(sq1) 
			da2  <- da1 * (pnorm((sq1*(1-theta) -qnorm(1-alpha2/2)*sqrt(2-2*theta))/sqrt(1-theta^2))-
                pnorm(( lam1 -qnorm(1-alpha2/2)*sqrt(2-2*theta)/2 -sq1*theta )/sqrt(1-theta^2)))
			Interseccio <- 2*(prob1  <- sum(da2*diff(sq1)[1])/(alpha2/2))
		})
		Inter 	<- mean(fas2)
		alpha2T <- (2*lambda1 -Inter)*alpha2
	}
	else{
		alpha1A <- lambda1
		alpha2T <- 0
	}
		
     ## finish initial iteration 
 	 Z             <- flsaFast(A, L=1, lambda1, lambda2, corsOmega = corsOmega, P=P, 
							 ind12=ind12, nsubset = nsubset, automLambdas = automLambdas,  penalize.diagonal = TRUE,
							 sigmaEstimate = sigmaEstimate, temporalFolders = FALSE, notOnlyLambda2 = notOnlyLambda2, 
							 alpha1A=alpha1A, isnet = FALSE)

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
           thetae[[k]]  <- theta[[k]] %*% (C12[[k]]*(N-1)/N  + (Z[[k]] - U[[k]]))
           A[[k]]      <- thetae[[k]] + U[[k]]
        }

 	    Z             <- flsaFast(A, L=1, lambda1, lambda2, corsOmega = corsOmega, P=P, 
							 ind12=ind12, nsubset = nsubset, automLambdas = automLambdas,  penalize.diagonal = TRUE,
							 sigmaEstimate = sigmaEstimate, temporalFolders = FALSE, notOnlyLambda2 = notOnlyLambda2, 
							 alpha1A=alpha1A, isnet = FALSE)

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
	obj$alpha2T  		<- alpha2T
		
    return(obj)
}


