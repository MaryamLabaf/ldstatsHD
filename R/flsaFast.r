######### thresholding step in joint estimation algorithm ###########
flsaFast <- function (A, L, lam1, lam2, corsOmega = NULL, P, ind12, 
                      nsubset = 10000, automLambdas = TRUE, penalize.diagonal = FALSE, 
                      sigmaEstimate = "CRmad", temporalFolders, notOnlyLambda2 = TRUE)
{
        ## lambda2 expression
        P		    <- dim(A[[1]])[1]
        q		    <- dim(A[[1]])[2]
        
        if(P == q) P2 			<- P * (P-1)/2
        else 	   P2 			<- P * q
        varsOmega 	<- sqrt(2 - 2 * corsOmega)
        AUX       	<- ((A[[1]]) - (A[[2]]))/varsOmega
                
        if(temporalFolders){
     		save(A, file = "Atemp.Rdata"); 
     		save(AUX, file = "Auxtemp.Rdata"); 
     		rm(A, AUX); gc()
   		 }
   		 
        if(automLambdas){
         if(is.null(nsubset)) id1 <- 1:P2
    	 else                 id1 <- sample(1:P2, min(P2,nsubset))
		 whichs    <- ind12[id1,]
          
         if(temporalFolders){
     		load("Auxtemp.Rdata"); 
     	 }
   		 if(sigmaEstimate == "mad")
             sigma     <- mean(abs(AUX[whichs]))/sqrt(2/pi)
         if(sigmaEstimate == "IQR")
             sigma <-  (quantile(AUX[whichs],0.75) - quantile(AUX[whichs],0.25))/1.349
         if(sigmaEstimate == "CRmad")
             sigma <-  s_Qn(AUX[whichs])
             
         lam2e     <- qnorm(1-lam2) * sigma
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
        rm(Y1,X1,S1,S2,Y2,X2,X3,Y3,S3,A);gc()

        ## lambda1 expression
        if(automLambdas&notOnlyLambda2){
          
          ddd     <- c(lowerTri(a1[whichs]),lowerTri(a2[whichs]))
          if(sigmaEstimate=="mad")
             sigma     <- mean(abs(ddd))/sqrt(2/pi)
          if(sigmaEstimate=="IQR")
             sigma <-  (quantile(ddd,0.75)-quantile(ddd,0.25))/1.349
          if(sigmaEstimate=="CRmad")
             sigma <-  s_Qn(ddd)
          rm(ddd); gc()
 
          lam1e       <- qnorm(1-lam1) * sigma
          lam1        <- lam1e
        }

        ## soft thresholding
        X = softA(a = a1, lam = lam1/L, penalize.diagonal = penalize.diagonal)
        rm(a1); gc()
        Y = softA(a = a2, lam = lam1/L, penalize.diagonal = penalize.diagonal)
        rm(a2); gc()
        
        ## temporal folders removal
        if(temporalFolders){
          file.remove("Auxtemp.Rdata")
          file.remove("Atemp.Rdata")
        }
    	return(list(Matrix(X), Matrix(Y)))
}

