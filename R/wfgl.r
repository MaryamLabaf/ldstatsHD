############# Joint partial correlation estimator ###########
wfgl <- function(D1, D2, lambda1, lambda2, paired = TRUE, automLambdas = TRUE, 
	 			 sigmaEstimate = "CRmad", pairedEst = "Reg-based-sim", maxiter = 30, 
	 			 tol = 1e-05, nsubset = 10000, weights = c(1,1), rho = 1, rho.increment = 1, 
	 			 triangleCorrection = TRUE, alphaTri = 0.01, temporalFolders = FALSE, 
	 			 notOnlyLambda2 = TRUE, roundDec = 16, burn = 0, lambda1B = NULL, lambda2B = NULL)		
{    
    
    ## checks
    cj <- controlwfgl(D1 = D1, D2 = D2, lambda1 = lambda1, lambda2 = lambda2, paired = paired, 
    			automLambdas = automLambdas, sigmaEstimate = sigmaEstimate, pairedEst = pairedEst, 
    			maxiter = maxiter, tol = tol, nsubset = nsubset, weights = weights, rho = rho,
	 			rho.increment = rho.increment, triangleCorrection = triangleCorrection, 
	 			alphaTri = alphaTri, temporalFolders = temporalFolders, burn = burn)		
	if(automLambdas) sigmaEstimate 	<- cj[1]
	if(paired) 		 pairedEst 		<- cj[2]
	 			 
    ## Initialisation
    lengthPath	  <- length(lambda1)
    D1            <- scale(D1)
    D2            <- scale(D2)
    N             <- dim(D1)[1]
    S			  <- list()
    S[[1]]        <- round(cov(D1)*(N - 1)/N,roundDec)
    S[[2]]        <- round(cov(D2)*(N - 1)/N,roundDec)
 	theta0        <- list()
    P             <- dim(S[[1]])[2]
    K             <- 2
    n             <- weights
 		       
	## theta0
	for(k in 1:K){
		aux 		<- S[[k]]
		edecomp     <- eigen(aux)
		rm(aux); 
		D           <- edecomp$values
		V           <- edecomp$vectors
		rm(edecomp); 
		DK          <- n[k]/(2 * rho) * (-D + sqrt(D^2 + 4 * rho/n[k]))
		AA          <- t(apply(V,1,function(x) x*DK))
		theta0[[k]] <- AA %*% t(V)
		rm(V); rm(DK); rm(AA); 
    }
	## corEst2
	if(paired){
		c12 		<- apply(as.matrix(1:P),1,function(d) theta0[[1]][d,]/(theta0[[1]][d,d]))
		c22 		<- apply(as.matrix(1:P),1,function(d) theta0[[2]][d,]/(theta0[[2]][d,d]))
		cors 	    <- -apply(as.matrix(1:P),1,function(j) cor(D1[,j] + D1[,-j] %*% c12[-j,j], D2[,j] + D2[,-j]%*%c22[-j,j]))
		pcor1 		<- cov2cor(theta0[[1]])
		pcor2 		<- cov2cor(theta0[[2]])
		if(pairedEst == "Reg-based-sim")
		{
			corEst2 	<- (cors%*%t(cors)) * (1/((1 - pcor1^2) * (1 - pcor2^2)))
			diag(corEst2) <- cors^2
		}
		if(pairedEst == "Reg-based")
		{
			W1 			<- array(rep(cors,P), dim = c(P, P))
			corEst2 	<- (cors%*%t(cors) + (pcor1 * pcor2 * (W1^2 + t(W1)^2) / 2)) * 
						   (1/((1 - pcor1^2) * (1 - pcor2^2)))
			diag(corEst2) <- cors^2
			rm(W1)
		}
		corEst2 <- pmin(corEst2,0.95)
		rm(cors,c12,c22); 
	}
	else  corEst2 <- array(0.5, dim=c(P,P))

	if(temporalFolders){
 	  save(theta0, file = "theta0temp.Rdata")
	  rm(theta0); 
	}
	
	## ind12
	if(automLambdas)  ind12 <- do.call(cbind,lowerTriInd(P))
	else              ind12 <- NULL

    if(length(lambda1) > 1){
  	 objaux <- list()
  	 for(lamb1i in 1:lengthPath){
  	 	if(temporalFolders)    load("theta0temp.Rdata")
  	 	
       objaux[[lamb1i]] <- wfglaux(D1 = D1, D2 = D2, lambda1 = lambda1[lamb1i], lambda2 = lambda2, paired = paired, 
								 automLambdas = automLambdas, sigmaEstimate = sigmaEstimate, pairedEst = pairedEst,
								 maxiter = maxiter, tol = tol, nsubset = nsubset, weights = weights, rho = rho, 
								 rho.increment = rho.increment, triangleCorrection = triangleCorrection, 
								 alphaTri = alphaTri, temporalFolders = temporalFolders, theta0 = theta0, 
								 corEst2 = corEst2, ind12 = ind12, S = S, notOnlyLambda2 = notOnlyLambda2, burn=burn,
	 							 lambda1B = lambda1B, lambda2B = lambda2B)	
	 }
	 obj 	    <- list()
	 obj$path   <- lapply(objaux, function(x) x$path)
	 obj$omega  <- lapply(objaux, function(x) x$omega)
	 obj$triangleCorrection <- lapply(objaux, function(x) x$triangleCorrection)
	 obj$weakTriangEdges <- lapply(objaux, function(x) x$weakTriangEdges)
	 obj$weakTriangEdgesPval <- lapply(objaux, function(x) x$weakTriangEdgesPval)
	 obj$iters <- lapply(objaux, function(x) x$iters)
	 obj$triangles1 <- lapply(objaux, function(x) x$triangles1)
	 obj$triangles2 <- lapply(objaux, function(x) x$triangles2)
	 obj$corEst2 	<- corEst2	
	 obj$automLambdas	<- automLambdas
	 rm(objaux)

	}
	else{
  	 	if(temporalFolders)    load("theta0temp.Rdata")
		obj <- wfglaux(D1 = D1, D2 = D2, lambda1 = lambda1, lambda2 = lambda2, paired = paired, 
					   automLambdas = automLambdas, sigmaEstimate = sigmaEstimate, pairedEst = pairedEst,
					   maxiter = maxiter, tol = tol, nsubset = nsubset, weights = weights, rho = rho, 
					   rho.increment = rho.increment, triangleCorrection = triangleCorrection, 
					   alphaTri = alphaTri, temporalFolders = temporalFolders, theta0 = theta0, 
 					   corEst2 = corEst2, ind12 = ind12, S = S, notOnlyLambda2 = notOnlyLambda2, burn=burn,
 					   lambda1B = lambda1B, lambda2B = lambda2B)	
	}
	
	obj$paired  		<- paired
	obj$sigmaEstimate  	<- sigmaEstimate
 	obj$pairedEst  		<- pairedEst
    obj$lambda1  		<- lambda1
	obj$lambda2  		<- lambda2
    obj$automLambdas	<- automLambdas

    class(obj) 			<- "wfgl"
    
   	if(temporalFolders){
	 file.remove("ind12temp.Rdata")
	 file.remove("corEst2temp.Rdata")
	 file.remove("stemp.Rdata")
	 file.remove("theta0temp.Rdata")
	 file.remove("d1d2temp.Rdata")
	}

    return(obj)
}

             
#######
wfglaux <- function(D1, D2, lambda1, lambda2, paired = TRUE, automLambdas = TRUE, 
	 			 sigmaEstimate = "CRmad", pairedEst = "Reg-based-sim", maxiter = 30, 
	 			 tol = 1e-05, nsubset = 10000, weights = c(1,1), rho = 1, rho.increment = 1, 
	 			 triangleCorrection = TRUE, alphaTri = 0.01, temporalFolders = FALSE, 
	 			 theta0 = NULL, corEst2 = NULL, ind12 = NULL, S = NULL, notOnlyLambda2 = TRUE, 
	 			 burn = 0, lambda1B = NULL, lambda2B = NULL)		
{	
	
	theta         <- theta0 
    if(temporalFolders){
     save(D1,D2,file="d1d2temp.Rdata")
     save(ind12,file="ind12temp.Rdata")
     save(corEst2,file="corEst2temp.Rdata")
     save(S,file = "stemp.Rdata")
     rm(D1,D2, S, theta0, corEst2, ind12);
    }
	W			  <- list()
	    
    P             <- dim(theta[[1]])[2]
    K             <- 2
    n             <- weights
    diff_value    <- numeric()
    sigmas 		  <- numeric()
    theta.prev    <- list()
    A			  <- list()

	############ initial iteration  
	iter          <- 0
	for (k in 1:K) {
    	W[[k]] 		    <- matrix(0, P, P)
		A[[k]]    	    <- theta[[k]] 
		theta.prev[[k]] <- matrix(0, P, P)
	}
	if(temporalFolders){
 			 save(theta, file= "thetatemp.Rdata")
 			 save(W, file= "wtemp.Rdata")
 			 save(theta.prev, file= "thetaPtemp.Rdata")
			 rm(theta, W, theta.prev);
			 load("ind12temp.Rdata")
			 load("corEst2temp.Rdata")
	}
	
	if(lambda2>0&automLambdas){
	    sq1  <- seq(-8,8,length.out=10000)
		alpha2 <- lambda2
		thetaseq <- seq(-0.3,0.99,length.out=1000)
		fas <- apply(as.matrix(thetaseq),1,function(theta){
         x2 <- sort(abs(sq1-qnorm(1-alpha2/2)*sqrt(2-2*theta)/2),index.return=TRUE)
		 x2$x[which(cumsum( (dnorm(sq1)*diff(sq1)[1] * pnorm((sq1*(1-theta) -qnorm(1-alpha2/2)*sqrt(2-2*theta))/sqrt(1-theta^2))/(alpha2/2))[x2$ix])>1-lambda1)[1]]
        })
		iddd 		<- pmax(1,round(punif(corEst2,-0.3,0.99)*1000))
		alpha1A 	<- matrix(fas[iddd],ncol=dim(corEst2)[1]) 
		
		thetaSeq2 <- sample(1:length(corEst2),1000)
 		fas2 <- apply(as.matrix(thetaSeq2),1,function(j){
 			theta <- corEst2[j]
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
	
	automLambdase 	<- automLambdas
	lambda1e 		<- lambda1
	lambda2e 		<- lambda2
	if(burn>0){
	 automLambdas <- FALSE
	 lambda1 	  <- lambda1B
	 lambda2 	  <- lambda2B
	} 	
	Z             <- flsaFast(A, rho, lambda1, lambda2, corsOmega = corEst2, P=P, 
							 ind12, nsubset = nsubset, automLambdas = automLambdas, sigmaEstimate = sigmaEstimate,
							 temporalFolders = temporalFolders, notOnlyLambda2 = notOnlyLambda2, alpha1A = alpha1A,
							 isnet = TRUE)
    if(temporalFolders){
				save(Z,file="ztemp.Rdata"); 
				load("thetatemp.Rdata")
				load("Wtemp.Rdata")
	}

	## Updating step
	for (k in 1:K) {
		W[[k]]    <- W[[k]] + (theta[[k]] - Z[[k]])
	}
	if(temporalFolders){
		save(W, file="wtemp.Rdata"); 
		load("thetaPtemp.Rdata");
		rm(Z, W); 
	}
		
	## Convergence check
	iter = iter + 1
	diff_value[iter] = 0
	for (k in 1:K) {
		diff_value[iter] = diff_value[iter] + sum(abs(theta[[k]] - theta.prev[[k]]))/sum(abs(theta.prev[[k]]))
	}
			if(temporalFolders){
			 rm(theta,theta.prev)
			}
			rho = rho * rho.increment

   
    ############# recursive algorithm
     pb    <- txtProgressBar(min = 0, max = maxiter, style = 3)
     while (iter <= maxiter && diff_value[iter] > tol) {    
		setTxtProgressBar(pb, iter)    
		if(temporalFolders){
		 load("thetatemp.Rdata")
		}
		theta.prev  <- theta
		if(temporalFolders){
			save(theta.prev, file = "thetaPtemp.Rdata"); 
			rm(theta.prev);
		 }
	
		## Maximisation step
		for (k in 1:K) {
		  if(temporalFolders){
		   load("stemp.Rdata")
		   load("Wtemp.Rdata")
		   load("Ztemp.Rdata")
		  }
		  aux <- S[[k]] - rho * Z[[k]]/n[k] + rho * W[[k]]/n[k]
		  if(temporalFolders){
		   rm(S,Z,W);gc()
		  }
		  edecomp     <- eigen(aux)
		  rm(aux);
		  D           <- edecomp$values
		  V           <- edecomp$vectors
		  rm(edecomp);
		  DK          <- n[k]/(2 * rho) * (-D + sqrt(D^2 + 4 * rho/n[k]))
		  AA          <- t(apply(V,1,function(x) x*DK))
		  theta[[k]]  <- AA %*% t(V)
		  rm(V); rm(DK); rm(AA)
		}
		if(temporalFolders){
			save(theta, file = "thetatemp.Rdata"); 
			rm(theta);gc()
		 }
		
 		 if(temporalFolders){
			 load("thetatemp.Rdata")
			 load("Wtemp.Rdata")
		 }
		
		 ## Thresholding step        
		 for (k in 1:K) {
			A[[k]]    <- theta[[k]] + W[[k]]
		 }
		 if(temporalFolders){
		    rm(theta,W);gc()
 		    load("ind12temp.Rdata")
			load("corEst2temp.Rdata")
			
		 }
		if(iter>burn){
		 automLambdas <- automLambdase
	     lambda1 	  <- lambda1e
	 	 lambda2 	  <- lambda2e
		} 
	    Z             <- flsaFast(A, rho, lambda1, lambda2, corsOmega = corEst2, P=P, 
							 ind12, nsubset = nsubset, automLambdas = automLambdas, sigmaEstimate = sigmaEstimate,
							 temporalFolders = temporalFolders, notOnlyLambda2 = notOnlyLambda2, alpha1A=alpha1A,isnet=TRUE)	
		 if(temporalFolders){
				rm(corEst2, ind12); gc()
				save(Z,file="ztemp.Rdata"); 
				load("thetatemp.Rdata")
				load("Wtemp.Rdata")
		 }
			## Updating step
			for (k in 1:K) {
				W[[k]]    <- W[[k]] + (theta[[k]] - Z[[k]])
			}
			 if(temporalFolders){
				save(W,file="wtemp.Rdata"); 
				load("thetaPtemp.Rdata");
				rm(Z,W);gc()
			}
		
			## Convergence check
			iter = iter + 1
			diff_value[iter] = 0
			for (k in 1:K) {
				diff_value[iter] = diff_value[iter] + sum(abs(theta[[k]] - theta.prev[[k]]))/sum(abs(theta.prev[[k]]))
			}
			if(temporalFolders){
			 rm(theta,theta.prev)
			}
			rho = rho * rho.increment

	}
	close(pb)
	if(temporalFolders){
		  load("ztemp.Rdata");
		  file.remove("wtemp.Rdata")
		  file.remove("thetatemp.Rdata")
		  file.remove("thetaPtemp.Rdata")
	 }
	## triangle correction
	A1                  <- as.matrix(Z[[1]])!=0 * 1 
	A2                  <- as.matrix(Z[[2]])!=0 * 1 
	diag(A1)            <- 0
	diag(A2)            <- 0
	if(temporalFolders){
	 rm(Z); gc()
	 load("d1d2temp.Rdata")
	}

	if(triangleCorrection & !is.null(D1) & !is.null(D2) ){
	  triang              <- pvaluesTriangleWeakestEdge(A1, A2, alpha = alphaTri, D1, D2)
	  graphstructure      <- list(Matrix(triang$A1.new), Matrix(triang$A2.new)) 
	  weakTriangEdges     <- triang$edgesToConsider
	  weakTriangEdgesPval <- triang$pvaluesToConsider
	  triangles1 		  <- triang$TRI1r 
	  triangles2 		  <- triang$TRI2r
	}
	else{
	 graphstructure      <- list(Matrix(A1),Matrix(A2))
	 weakTriangEdges     <- NULL
	 weakTriangEdgesPval <- NULL
	  triangles1 		 <- NULL
	  triangles2 		 <- NULL
	}
   
	if(temporalFolders){
	 load("ztemp.Rdata");
	 file.remove("ztemp.Rdata")
	 load("ind12temp.Rdata")
	 load("corEst2temp.Rdata")
	 load("stemp.Rdata")
	}

    obj <- list(path = graphstructure, omega = Z, triangleCorrection = triangleCorrection,
               weakTriangEdges = weakTriangEdges, weakTriangEdgesPval = weakTriangEdgesPval, 
               diff_value = diff_value, iters = iter, corEst = corEst2, triangles1 = triangles1,
               triangles2 = triangles2)

	obj$paired  		<- paired
	obj$sigmaEstimate  	<- sigmaEstimate
	obj$pairedEst  		<- pairedEst
	obj$lambda1  		<- lambda1
	obj$lambda2  		<- lambda2
	obj$alpha2T  		<- alpha2T
	obj$alpha1A  		<- alpha1A
	
	return(obj)

}


