
########################## Partial correlation and data simulator ##############
pcorSimulatorJoint  <- function(nobs, nclusters, nnodesxcluster, pattern = "hubs", 
								diffType = "cluster", dataDepend = "ind", low.strength = 0.5, 
								sup.strength = 0.9, pdiff = 0, nhubs = 5, degree.hubs = 20,  
								nOtherEdges = 30, alpha = 2.3, plus = 0, prob = 0.05, 
								perturb.clust = 0, mu = 0, diagCCtype = "dicot", 
								diagNZ.strength = .5, mixProb = 0.5, probSign = 0.5,  
								exactZeroTh = 0.05, seed = sample(10000,nclusters+2))
{

       ## Checks
       cJPS <- controlsPcorSimulatorJoint(nobs, nclusters, nnodesxcluster, pattern, 
                      diffType, dataDepend, low.strength, sup.strength, pdiff, nhubs, 
                      degree.hubs, nOtherEdges, alpha, plus, prob, perturb.clust, mu, 
                      diagCCtype, diagNZ.strength, mixProb, probSign, seed, exactZeroTh)
                      
        ## two similar precision matrices
        PCOR   <- try(pcorSimulatorTwo(nobs = cJPS$nobs, nclusters = cJPS$nclusters, 
        			  nnodesxcluster = cJPS$nnodesxcluster, pattern = cJPS$pattern, 
        			  diffType = cJPS$diffType, low.strength = cJPS$low.strength, 
        			  sup.strength = cJPS$sup.strength, nhubs = cJPS$nhubs, 
        			  degree.hubs = cJPS$degree.hubs, nOtherEdges = cJPS$nOtherEdges, 
        			  alpha = cJPS$alpha, plus = cJPS$plus, seed = cJPS$seed, prob = cJPS$prob, 
        			  perturb.clust = cJPS$perturb.clust, mu = cJPS$mu, pdiff = cJPS$pdiff,  
        			  mixProb = cJPS$mixProb, probSign = cJPS$probSign, exactZeroTh = cJPS$exactZeroTh))
       
        M1     <- PCOR$omega0
        M2     <- PCOR$omega1
        P      <- dim(M1)[1]
        N      <- nobs
        
        ## Relevant differential edges
        PCOR1 <- M1
        if(pdiff == 0) M2 <- M1
        PCOR2 <- M2
        DIFlargThr <- abs(cov2cor(solve(PCOR1)) - cov2cor(solve(PCOR2))) > cJPS$exactZeroTh
        diffs      <- (PCOR1-PCOR2 !=0) * DIFlargThr
        
       
       if(cJPS$dataDepend == "ind"){
            CC                  <- 1
            delta               <- 0 
            jointPCOR0          <- as.matrix(bdiag(PCOR1,PCOR2))
        	omega0 				<- jointPCOR0[1:P,1:P]
		    omega1 				<- jointPCOR0[1:P+P,1:P+P]

        }
        else
        { 
			## paired structure 
			if(cJPS$diagCCtype == "dicot") 
			  delta  <- sample(c(rep(0,P/2),rep(diagNZ.strength,P/2)))
			if(cJPS$diagCCtype == "beta13") 
			  delta <- rbeta(P,1,3) 
				
			## diagonal in the partial correlation 
			if(cJPS$dataDepend == "diagOmega"){
				CC                  <- 1
				DAT                 <- 0 
				jointPCOR0          <- as.matrix(bdiag(PCOR1,PCOR2))
    			jointPCOR0[1:P,P+1:P]   <- -(diag(delta))
				jointPCOR0[1:P+P,1:P]   <- -(diag(delta))
				eigval 					<- 	eigen(jointPCOR0)$val
				if(class(eigval)=="complex") eigval <- Re(eigval)

            	while(max(eigval)/min(eigval) > 2*P| min(eigval) <0){
				  	diag(jointPCOR0)        <- diag(jointPCOR0) + 0.1
  					eigval 					<- 	eigen(jointPCOR0)$val
					if(class(eigval)=="complex") eigval <- Re(eigval)
				}
 				omega0 		<- jointPCOR0[1:P,1:P]
		        omega1 		<- jointPCOR0[1:P+P,1:P+P]
 			 }
 
			if(cJPS$dataDepend == "mult"){
				CC                    <- 1
				DAT                   <- 0 
				jointPCOR0            <- as.matrix(bdiag(PCOR1,PCOR2))
  			    COVj                <- cov2cor(solve(jointPCOR0))            
             	A1                  <- cov2cor(COVj[1:P,1:P])
				EIG                 <- eigen(A1)
				A1                  <- EIG$vectors%*% diag(sqrt(EIG$val))%*%t(EIG$vectors)
				A1                  <- apply(A1,1,as.numeric)
				A2                  <- cov2cor(COVj[1:P+P,1:P+P])
				EIG                 <- eigen(A2)
				A2                  <- EIG$vectors%*%diag(sqrt(EIG$val))%*%t(EIG$vectors)
				A2                  <- apply(A2,1,as.numeric)
				A                   <- A2%*%A1
				COVj[1:P,P+1:P]     <- diag(sqrt(delta))%*%A%*%diag(sqrt(delta))
				COVj[1:P+P,1:P]     <- diag(sqrt(delta))%*%A%*%diag(sqrt(delta))
				diag(COVj)          <- diag(COVj)
				  
				eigval 				<- eigen(COVj)$val
				if(class(eigval)=="complex") eigval <- Re(eigval)
				
				while(max(eigval)/min(eigval) > 2*P| min(eigval) <0)
			    {
				  diag(jointPCOR0)    <- diag(jointPCOR0) + 0.1
				  COVj                <- cov2cor(solve(jointPCOR0))            
             	  A1                  <- cov2cor(COVj[1:P,1:P])
				  EIG                 <- eigen(A1)
				  eigval 			  <- EIG$val
				  eigvec 			  <- EIG$vec
				  if(class(eigval)=="complex"){
				   eigval <- Re(eigval)
				   eigvec <- Re(eigvec)
				  } 
				  A1                  <- eigvec%*% diag(sqrt(eigval))%*%t(eigvec)
				  A1                  <- apply(A1,1,as.numeric)
				  A2                  <- cov2cor(COVj[1:P+P,1:P+P])
				  
				  EIG                 <- eigen(A2)
				  eigval 			  <- EIG$val
				  eigvec 			  <- EIG$vec
				  if(class(eigval)=="complex"){
				   eigval <- Re(eigval)
				   eigvec <- Re(eigvec)
				  } 
				  A2                  <- eigvec%*% diag(sqrt(eigval))%*%t(eigvec)
				  A2                  <- apply(A2,1,as.numeric)
				  A                   <- A2%*%A1
				  COVj[1:P,P+1:P]     <- diag(sqrt(delta))%*%A%*%diag(sqrt(delta))
				  COVj[1:P+P,1:P]     <- diag(sqrt(delta))%*%A%*%diag(sqrt(delta))
				  diag(COVj)          <- diag(COVj)
	   			  eigval 				<- eigen(COVj)$val
				  if(class(eigval)=="complex") eigval <- Re(eigval)					
				}
				omega0 				<- jointPCOR0[1:P,1:P]
		        omega1 				<- jointPCOR0[1:P+P,1:P+P]
				jointPCOR0 			<- solve(COVj)
			}
        
			if(cJPS$dataDepend == "add"){
				CC                  <- 1
				DAT                 <- 0 
			    jointPCOR0           <- as.matrix(bdiag(PCOR1,PCOR2))
				COVj                <- cov2cor(solve(jointPCOR0))            
				COVj[1:P,P+1:P]     <- diag(sqrt(delta))%*%COVj[1:P,1:P]%*%diag(sqrt(delta))
				COVj[1:P+P,1:P]     <- diag(sqrt(delta))%*%COVj[1:P,1:P]%*%diag(sqrt(delta))
				
				eigval 				<- eigen(COVj)$val
				if(class(eigval)=="complex") eigval <- Re(eigval)
				
				while(max(eigval)/min(eigval) > 2*P| min(eigval) <0)
			    {
				  diag(jointPCOR0)    <- diag(jointPCOR0) + 0.1
				  COVj                <- cov2cor(solve(jointPCOR0))            
				  COVj[1:P,P+1:P]     <- diag(sqrt(delta))%*%COVj[1:P,1:P]%*%diag(sqrt(delta))
				  COVj[1:P+P,1:P]     <- diag(sqrt(delta))%*%COVj[1:P,1:P]%*%diag(sqrt(delta))
 				  eigval 			  <- eigen(COVj)$val
				  if(class(eigval)=="complex") eigval <- Re(eigval)
				}

        	    omega0 				<- jointPCOR0[1:P,1:P]
		        omega1 				<- jointPCOR0[1:P+P,1:P+P]
		    }
		}

        COVj <- cov2cor(solve(jointPCOR0))
		DAT  <- try(mvrnorm(N,rep(0,P*2),COVj))
        D1 	 <- (DAT[,1:P])
        D2 	 <- (DAT[,1:P+P])

		path1        	<- ((omega0)!=0)*1
	    diag(path1)  	<- 0

		path2        	<- ((omega1)!=0)*1
	    diag(path2)  	<- 0

       obj <- list(D1 = D1, D2 = D2, omega1 = omega0, omega2 = omega1, 
                    P = P, diffs = diffs, delta = delta, covJ = COVj, dataDepend = cJPS$dataDepend, 
                    diagCCtype = cJPS$diagCCtype, pattern = cJPS$pattern, path1 = path1, path2 = path2)
        class(obj) <- "pcorSimJoint"          
        return(obj)
}




