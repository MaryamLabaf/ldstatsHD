thesholdSelection <- function(D1,D2, useq, deltaA=3, deltaB=10, nite= 500, excAdj=FALSE, alpha=0.05, paired= TRUE)
{				
		        if(paired) thetaKnown <- 0
		        else thetaKnown <- NULL
		        
		        Test  	<- ztransfCorrDiff(D1, D2, thetaKnown = thetaKnown )$Ts
			    pval 	<- c(2-2*pnorm(abs(Test)))
			    phis 	<- 1-pi0est(pval,0.30)$pi0		 
				M 		<- length(Test)
				N		<- dim(D1)[1]
				if(phis>0){
				delta <- numeric()
				SS   <- numeric()
				for (i in 1:nite){
				  SS 	<- c(SS,round(M*phis))
				  delta <- c(delta,abs(rgamma(SS[i],deltaA,deltaB)))
				}
				s2 <- c(0,cumsum(SS))
				power <- array(0,dim=c(length(useq),nite))
   				for(kk in 1:length(useq)){
   				    u <- useq[kk]
  				    #print(u)
  				    phi0  <- 2*(1-pnorm(u))
				    dt    <- delta*sqrt(N-3)
					D <- dt
					dn1 <- dnorm(u-D)
					pn1 <- pnorm(-(u-D))
					dn2 <- dnorm(-u-D)
					pn2 <- pnorm(-u-D)			    
				    phi1 <- pn1 + pn2
				    
				    if(excAdj){
				     mu01 	<- u^2 + 1 -u *dnorm(u)/pnorm(-u)			    	
			    	 var01 	<- 3 + u^4 +6*u^2 - (5*u +u^3) *dnorm(u)/pnorm(-u) - mu01^2			
				    
 				     varr01 <- M*phi0*((1-phi0)*mu01^2 +var01)
				     varr11 <- (M-SS)*phi0*((1-phi0)*mu01^2 +var01)
					(EXP10 <- (u-D) * dn1/(pn1+pn2) + 1 + D^2 + 2*D*(dn1-dn2)/(pn1+pn2) - (-u-D)*dn2/(pn2+pn1))
					
					(EXP11 <- EXP10 +u^2 - 2*u*(dn1+dn2)/(pn1+pn2) - 2*u*D*(pn1-pn2)/(pn1+pn2))
				     B <- (u-D)^2*dn1  +2*dn1
					 B2 <-  (u+D)^2*dn2  +2*dn2
					
					 A1 <- (4*u*(u-D)^2*dn1 + 8*u*dn1 + 12*D*u*(u-D)*dn1 +12*D*u*pn1 +12*u*D^2*dn1 + 4*u*D^3*pn1)/(pn1+pn2)
					 A2 <- (4*u*(-u-D)^2*dn2 + 8*u*dn2 + 12*D*u*(-u-D)*dn2 -12*D*u*pn2 +12*u*D^2*dn2 - 4*u*D^3*pn2)/(pn1+pn2)
	 				
					 B <- (u-D)^2*dn1  +2*dn1
				 	 B2 <-  -(-u-D)^2*dn2  -2*dn2
				 	 (EX210 <- ( D^4 + 4*D^3 *(dn1-dn2)/(pn1+pn2) + 6*D^2*(pn1 + pn2 + (((u-D)*dn1) -(-u-D)*dn2))/(pn1+pn2) + 4*D*(B+B2)/(pn1+pn2) + 
					 ((u-D)^3*dn1 + 3*(u-D)*dn1 + 3*pn1)/((pn1+pn2)) +( -(-u-D)^3*dn2 + 3*(u+D)*dn2 + 3*pn2)/((pn1+pn2))))
					  
					(Var11 <- EX210 -A1 -A2 +6*u^2*EXP10 -4*u^3 *(dn1+dn2)/(pn1+pn2) - 4*u^3*D *(pn1-pn2)/(pn1+pn2)  + u^4 -EXP11^2)


		 		     asd1 <- apply(as.matrix(2:(nite+1)),1,function(j){
				     	Num   <- sum(phi1[s2[j-1]:s2[j]]*EXP11[s2[j-1]:s2[j]]) - SS[j-1] * phi0 *mu01 - qnorm(1-alpha/2) * sqrt(varr01)
				     	Varr1 <- sum(phi1[s2[j-1]:s2[j]]*((1-phi1[s2[j-1]:s2[j]])*EXP11[s2[j-1]:s2[j]]^2 +Var11[s2[j-1]:s2[j]]))
						Den   <- sqrt(pmax( Varr1 + varr11[j-1],1))				  
  			        	(Num/Den)	
  			        })
				    
				   }
				    else{
					 mu00   <- 1 + u*dnorm(u)/(1-pnorm(u))
				     var00  <- (u^3+3*u) * dnorm(u)/(1-pnorm(u))+3 -mu00^2
		            
				     varr00 <- M*phi0*((1-phi0)*mu00^2 +var00)
				     varr10 <- (M-SS)*phi0*((1-phi0)*mu00^2 +var00)
				    
				     EXP10 <-(u-D) * dn1/(pn1+pn2) + 1 + D^2 + 2*D*(dn1-dn2)/(pn1+pn2) - (-u-D)*dn2/(pn2+pn1)
					 (Var10 <- ( D^4 + D^3 *(dn1-dn2)/(pn1+pn2) + 6*D^2 + D^2 * (u*dn2 + u*dn1)/(pn1+pn2) + D*(u^2+5)*(dn1-dn2)/(pn1+pn2) + 
					 (u^3+3*u)*(dn1+dn2)/(pn1+pn2) +3 - EXP10^2))

				     asd1 <- apply(as.matrix(2:(nite+1)),1,function(j){
				     	Num   <- sum(phi1[s2[j-1]:s2[j]]*EXP10[s2[j-1]:s2[j]]) - SS[j-1] * phi0 *mu00 - qnorm(1-alpha/2) * sqrt(varr00)
				     	Varr1 <- sum(phi1[s2[j-1]:s2[j]]*((1-phi1[s2[j-1]:s2[j]])*EXP10[s2[j-1]:s2[j]]^2 +Var10[s2[j-1]:s2[j]]))
						Den   <- sqrt(pmax( Varr1 + varr10[j-1],1))				  
  			        	(Num/Den)	
				    })
				   }
				   
				   power[kk,] <- as.numeric(asd1)
				 }
				   
				 bestThreshold <- useq[which.min(apply(apply(-power,2,rank),1,median))]				  
				}
				else 
				bestThreshold <- max(useq)
				return(bestThreshold)
}