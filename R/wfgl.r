############# Joint partial correlation estimator ###########
wfgl <- function(D1, D2, lambda1, lambda2, paired = TRUE, automLambdas = TRUE, 
	 			 sigmaEstimate = "CRmad", pairedEst = "Reg-based-sim", maxiter = 30, 
	 			 tol = 1e-05, nsubset = 10000, weights = c(1,1), rho = 1, rho.increment = 1, 
	 			 triangleCorrection = TRUE, alphaTri = 0.01, temporalFolders = FALSE, 
	 			 notOnlyLambda2 = TRUE)
{    
    
    ## checks
    cj <- controlwfgl(D1 = D1, D2 = D2, lambda1 = lambda1, lambda2 = lambda2, paired = paired, 
    			automLambdas = automLambdas, sigmaEstimate = sigmaEstimate, pairedEst = pairedEst, 
    			maxiter = maxiter, tol = tol, nsubset = nsubset, weights = weights, rho = rho,
	 			rho.increment = rho.increment, triangleCorrection = triangleCorrection, 
	 			alphaTri = alphaTri, temporalFolders = temporalFolders)
	if(automLambdas) sigmaEstimate 	<- cj[1]
	if(paired) 		 pairedEst 		<- cj[2]
	 			 
    ## Initialisation
    lengthPath	  <- length(lambda1)
    D1            <- scale(D1)
    D2            <- scale(D2)
    N             <- dim(D1)[1]
    S			  <- list()
    S[[1]]        <- cov(D1)
    S[[2]]        <- cov(D2)
 	theta0        <- list()
    P             <- dim(S[[1]])[2]
    K             <- 2
    n             <- weights
 		       
	## theta0
	for(k in 1:K){
		aux 		<- S[[k]]
		edecomp     <- eigen(aux)
		rm(aux); gc()
		D           <- edecomp$values
		V           <- edecomp$vectors
		rm(edecomp); gc()
		DK          <- n[k]/(2 * rho) * (-D + sqrt(D^2 + 4 * rho/n[k]))
		AA          <- t(apply(V,1,function(x) x*DK))
		theta0[[k]] <- AA %*% t(V)
		rm(V); rm(DK); rm(AA); gc()
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
			corEst2 	<- (cors%*%t(cors)) * (1/(sqrt(1 - pcor1^2) * sqrt(1 - pcor2^2)))
			diag(corEst2) <- cors^2
		}
		if(pairedEst == "Reg-based")
		{
			W1 			<- array(rep(cors,P), dim = c(P, P))
			corEst2 	<- (cors%*%t(cors) + (pcor1 * pcor2 * (W1^2 + t(W1)^2) / 2)) * 
						   (1/(sqrt(1 - pcor1^2) * sqrt(1 - pcor2^2)))
			diag(corEst2) <- cors^2
			rm(W1)
		}
		rm(cors,c12,c22); gc()
	}
	else  corEst2 <- array(.5, dim=c(P,P))

	if(temporalFolders){
 	  save(theta0, file = "theta0temp.Rdata")
	  rm(theta0); gc()
	}
	
	## ind12
	if(automLambdas)  ind12 <- do.call(cbind,lowerTriInd(P))
	else              ind12 <- NULL

    if(length(lambda1) > 1){
  	 objaux <- list()
  	 for(lamb1i in 1:lengthPath){
  	 	if(temporalFolders)    load("theta0temp.Rdata")
  	 	
    	pal           <- max(10^5,P*(P-1)/2)	   
        if(automLambdas & notOnlyLambda2& pal*lambda1[lamb1i]*lambda2>5)   
          lambda2 <- findingLambda2(pal, lam1=lambda1[lamb1i], lam2=lambda2, tol = 0.001, maxiter = 30) #internal use only

       objaux[[lamb1i]] <- wfglaux(D1 = D1, D2 = D2, lambda1 = lambda1[lamb1i], lambda2 = lambda2, paired = paired, 
								 automLambdas = automLambdas, sigmaEstimate = sigmaEstimate, pairedEst = pairedEst,
								 maxiter = maxiter, tol = tol, nsubset = nsubset, weights = weights, rho = rho, 
								 rho.increment = rho.increment, triangleCorrection = triangleCorrection, 
								 alphaTri = alphaTri, temporalFolders = temporalFolders, theta0 = theta0, 
								 corEst2 = corEst2, ind12 = ind12, S = S, notOnlyLambda2 = notOnlyLambda2)	
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
	 rm(objaux)

	}
	else{
  	 	if(temporalFolders)    load("theta0temp.Rdata")
		pal           <- max(10^5,P*(P-1)/2)	   
        if(automLambdas & notOnlyLambda2& pal*lambda1*lambda2>5)   
          lambda2 <- findingLambda2(pal, lam1=lambda1, lam2=lambda2, tol = 0.001, maxiter = 30) #internal use only
		obj <- wfglaux(D1 = D1, D2 = D2, lambda1 = lambda1, lambda2 = lambda2, paired = paired, 
					   automLambdas = automLambdas, sigmaEstimate = sigmaEstimate, pairedEst = pairedEst,
					   maxiter = maxiter, tol = tol, nsubset = nsubset, weights = weights, rho = rho, 
					   rho.increment = rho.increment, triangleCorrection = triangleCorrection, 
					   alphaTri = alphaTri, temporalFolders = temporalFolders, theta0 = theta0, 
 					   corEst2 = corEst2, ind12 = ind12, S = S, notOnlyLambda2 = notOnlyLambda2)	
	}
	
	obj$paired  		<- paired
	obj$sigmaEstimate  	<- sigmaEstimate
 	obj$pairedEst  		<- pairedEst
    obj$lambda1  		<- lambda1
	obj$lambda2  		<- lambda2

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
	 			 theta0 = NULL, corEst2 = NULL, ind12 = NULL, S = NULL, notOnlyLambda2 = TRUE)	
{

	theta         <- theta0 
    if(temporalFolders){
     save(D1,D2,file="d1d2temp.Rdata")
     save(ind12,file="ind12temp.Rdata")
     save(corEst2,file="corEst2temp.Rdata")
     save(S,file = "stemp.Rdata")
     rm(D1,D2, S, theta0, corEst2, ind12);gc();
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
			 rm(theta, W, theta.prev); gc();
			 load("ind12temp.Rdata")
			 load("corEst2temp.Rdata")
	}
	Z             <- flsaFast(A, rho, lambda1, lambda2, corsOmega = corEst2, P=P, 
							 ind12, nsubset = nsubset, automLambdas = automLambdas, sigmaEstimate = sigmaEstimate,
							 temporalFolders = temporalFolders, notOnlyLambda2 = notOnlyLambda2)	
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
		rm(Z, W); gc()
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
			rm(theta.prev);gc()
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
		  rm(aux);gc()
		  D           <- edecomp$values
		  V           <- edecomp$vectors
		  rm(edecomp);gc()
		  DK          <- n[k]/(2 * rho) * (-D + sqrt(D^2 + 4 * rho/n[k]))
		  AA          <- t(apply(V,1,function(x) x*DK))
		  theta[[k]]  <- AA %*% t(V)
		  rm(V); rm(DK); rm(AA)
		  gc()
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
		 Z      <- flsaFast(A, rho, lambda1, lambda2, corsOmega = corEst2, P=P, 
							 ind12, nsubset = nsubset, automLambdas = automLambdas, sigmaEstimate = sigmaEstimate,
							 temporalFolders = temporalFolders, notOnlyLambda2 = notOnlyLambda2)	
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

	return(obj)

}


findingLambda2 <- function(n, lam1, lam2, tol, maxiter=100)
{
a1 <- rnorm(n,0,1)
a2 <- rnorm(n,0,1)
thr1    <- 1-lam2
thrm    <- 0.5
thrM    <- 1
ite <- 1
NOTYET  <- TRUE
 while(NOTYET|ite>maxiter){
		thr <- quantile(abs(a1-a2),thr1) 
        S1      <-  abs(a1-a2) <= thr
        S2      <-  (a1-a2) >  thr
        S3      <- -(a1-a2) >  thr
     	X1      <- (a1+a2)/2
        Y1      <- (a1+a2)/2
        X2      <- a1 - thr/2
        Y2      <- a2 + thr/2
        X3      <- a1 + thr/2
        Y3      <- a2 - thr/2
        a11      <- (S1 * X1 + S2 * X2 + S3 * X3)
        a22      <- (S1 * Y1 + S2 * Y2 + S3 * Y3)
        th2 <- quantile(c(a11,a22),1-lam1)
	    X1      <- a11 - th2
        X2      <- a11 + th2
	    Y1      <- a22 - th2
        Y2      <- a22 + th2
        S11      <- a11 > th2
        S12      <- a11 < -th2
        S21      <- a22 > th2
        S22      <- a22 < -th2
        X = S11*X1 + S12*X2
        Y = S21*Y1 + S22*Y2        
        (diffr <- mean(X-Y!=0)/mean(X!=0|Y!=0) - lam2 )
        if(abs(diffr) <tol) NOTYET <- FALSE
        else{
        if(diffr >0)
        {
          thrm <- thr1
          thr1 <- (thr1 + thrM)/2
          
        }
        if(diffr < 0)
        {
          thrM <- thr1
          thr1 <- (thr1 + thrm)/2

        }
        }
        (ite <- ite+1)
}
return(1-thr1) 
}      
