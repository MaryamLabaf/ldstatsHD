######### thresholding step in joint estimation algorithm ###########
flsaFast <- function (A, L, lam1, lam2, corsOmega = NULL, P, ind12, 
                      nsubset = 10000, automLambdas = TRUE, penalize.diagonal = FALSE, 
                      sigmaEstimate = "CRmad", temporalFolders, notOnlyLambda2 = TRUE, 
                      isnet = TRUE, alpha1A = NULL){
        ## lambda2 expression
        P		    <- dim(A[[1]])[1]
        q		    <- dim(A[[1]])[2]
        
        if(isnet) P2 			<- P * (P-1)/2
        else 	   P2 			<- P * q
        varsOmega 	<- sqrt(2 - 2 * corsOmega)        
        AUX       	<- ((A[[1]]) - (A[[2]]))/varsOmega
        if(automLambdas&isnet){
        		sigmae  <-  s_Qn(lowerTri(A[[1]]))
				sigmae2 <-  s_Qn(lowerTri(A[[2]]))
				sigma	<- (sigmae+sigmae2)/2
				sigmae3 <-  s_Qn(lowerTri(AUX))
		}
        if(automLambdas&!isnet){
        		sigmae  <-  s_Qn(A[[1]])
				sigmae2 <-  s_Qn(A[[2]])
				sigma	<- (sigmae+sigmae2)/2
				sigmae3 <-  s_Qn(AUX)
		}
        if(temporalFolders){
     		save(A, file = "Atemp.Rdata"); 
     		rm(A); gc()
   		 }
   		 
        if(automLambdas){
         lam2e     <- pmin(qnorm(1-lam2/2) * sigmae3,99990)
         lam2      <- lam2e/2
        }
        
        ## similarity thresholding
        S1      <-  abs(AUX) <= lam2/L*2
        S2      <-  AUX >  2 * lam2
        S3      <- -AUX >  2 * lam2
        if(temporalFolders){
            rm(AUX);gc()
     		load("Atemp.Rdata"); 
     	}
     	
     	X1      <- (A[[1]] + A[[2]])/2
        Y1      <- X1
        X2      <- A[[1]] - (lam2)/L * varsOmega
        Y2      <- A[[2]] + (lam2)/L * varsOmega
        X3      <- A[[1]] + (lam2)/L * varsOmega
        Y3      <- A[[2]] - (lam2)/L * varsOmega
        a1      <- (S1 * X1 + S2 * X2 + S3 * X3)
        a2      <- (S1 * Y1 + S2 * Y2 + S3 * Y3)		
        rm(Y1,X1,S2,Y2,X2,X3,Y3,S3,A)

        ## lambda1 expression
        if(automLambdas&notOnlyLambda2){
          
          lam1e      <- pmin(qnorm(1-lam1/2) * sigma * sqrt(2+2*corsOmega)/2,99990)#s_Qn(a1[S1]),99990) #
#            qnorm(1-lam1/2, 0, sd(as.numeric(a1)))
        }
        else{
         alpha1A <- alpha1A
         lam1e <- alpha1A
         sigma <- 1
        }

        ## soft thresholding
        X = softA(a = a1, lam = lam1e/L * S1  + alpha1A * sigma * (1-S1), penalize.diagonal = penalize.diagonal)
        rm(a1); 
        Y = softA(a = a2, lam = lam1e/L * S1  + alpha1A * sigma * (1-S1), penalize.diagonal = penalize.diagonal)
        rm(a2); 
        
        ## temporal folders removala
        if(temporalFolders){
          file.remove("Atemp.Rdata")
        }
        
    	return(list(Matrix(X), Matrix(Y)))
}


