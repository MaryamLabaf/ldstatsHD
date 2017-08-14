#####################
equalCorrelationsTestExc <- function(Test, theta = NULL, TestsP = NULL, EX = 2.3, excAdj = TRUE,
									 dependency = TRUE, conf.level = 0.95, nite = 200, MINint =2, 
									 MAXint=100, saddlePoint = FALSE, ...)
{
      pgpd2 	<- fExtremes::pgpd
      aT    	<- Test#rnorm(length(Test))#
      P2    	<- length(Test)
      weight1 	<- ifelse(excAdj, 1, 0)
      
      if(dependency)
      {

         if(length(EX)>1){
            PVALandCI  <- t(apply(as.matrix(EX),1,function(DA){
    		 	 Exc     	<- squareWeight(abs(aT), DA, weight1)
			     maxsR2   	<- apply(abs(TestsP), 1, function(x) squareWeight(abs(x), DA, weight1))
            	 mu 		<- mean(maxsR2)
			     sigma2 	<- var(maxsR2)
  		         pvalue  	<- 1- pgamma(Exc, mu^2/sigma2, mu/sigma2)
  		         ci  		<- Exc - qnorm(c(conf.level, 0.0001), mu, sqrt(sigma2))
	    	     maxsT  	<- Exc - mu
	    	     c(pvalue,ci,maxsT)
            }))
            pvalue <- PVALandCI[,1]
            ci     <- PVALandCI[,c(2,3)]
            maxsT  <- PVALandCI[,4]
         }

         else{
             DA 		<- EX
   	     	 u			<- DA
   	     	 Exc     	<- squareWeight(abs(aT), DA, weight1)
             maxsR2   	<- apply(abs(TestsP), 1, function(x) squareWeight(abs(x), DA, weight1))
             
			 mu 		<- mean(maxsR2)
			 sigma2 	<- var(maxsR2)
  		     pvalue  	<- 1- pgamma(Exc, mu^2/sigma2, mu/sigma2)
  		     
     	     ci  		<- t(as.matrix(Exc - qnorm(c(conf.level, 0.0001), mu, sqrt(sigma2))))
    	     maxsT  	<- Exc - mu
         }

      }
      else{
        if(length(EX)>1){
            PVALandCI  <- t(apply(as.matrix(EX), 1, function(DA){
				 Exc     	<- squareWeight(abs(aT), DA, weight1)
				 mcarlo 	<- apply(as.matrix(1:nite),1,function(j) squareWeight(abs(rnorm(P2,0,1)), DA, weight1))
				 mu0    	<- mean(mcarlo)
            	 sigma2 	<- var(mcarlo)
				 k	    	<- mu0^2/sigma2
				 b      	<- sigma2/mu0
				 Res2  		<- 1 - (pgamma(Exc, k, 1/b))
				 Res3  		<- Exc - qgamma(c(conf.level, 0.0001), k, 1/b)
    	         c(Res2, Res3, Exc - mu0)
            }))
            pvalue <- PVALandCI[,1]
            ci     <- PVALandCI[,2:3]
            maxsT  <- PVALandCI[,4]
         }	
		else{
            u			<- EX
            Exc     	<- squareWeight(abs(aT), u, weight1)
            NeR     	<- sum(u)
			numberof    <- MAXint-MINint+1
			M			<- length(aT)
            IIF <- ifelse(MINint==1,2,1)
			if(excAdj){
             Ex 		<- u^2 + 1 -u *dnorm(u)/pnorm(-u)
			 VarX 		<- 3 + u^4 +6*u^2 - (5*u +u^3) *dnorm(u)/pnorm(-u) - Ex^2	
			 mu     	<- M * 2*pnorm(-u) * Ex
  			 sigma2 	<- M * 2*pnorm(-u) *( (1-2*pnorm(-u)) * Ex^2 +VarX)			
		     pvalue	 	<- 1 - pnorm(Exc, mu, sqrt(sigma2))			 
   			 ci  		<- t(as.matrix(Exc - qnorm(c(conf.level, 0.0001), mu, sqrt(sigma2))))
    	     maxsT  	<- Exc - mu
    	    }
            else{
			  Ex 	<- 1 +u *dnorm(u)/pnorm(-u)
			  VarX 	<- (u^3 +3*u)* dnorm(u)/pnorm(-u) +3 - Ex^2		
			  mu     <- M * 2*pnorm(-u) * Ex
  			  sigma2 <- M * 2*pnorm(-u) *( (1-2*pnorm(-u)) * Ex^2 +VarX)				 	
			  if(saddlePoint)
              {
                pvalue <- saddlePointFunction(Exc, NeR, M, u, ttmin = -1.45, ttmax = .4999, lengthTT = 5000, kmin = MINint, kmax = MAXint)
                pvalue <- 1- sum(pvalue,na.rm=TRUE)
              }  
			  else
			    pvalue <- 1 - pnorm(Exc, mu, sqrt(sigma2))			 
     	      ci  		<- t(as.matrix(Exc - qnorm(c(conf.level, 0.0001), mu, sqrt(sigma2))))
    	      maxsT  		<- Exc - mu
		   }
		}
		}	 			      
        return(list(testExc = maxsT, pvalue = pvalue, ci = ci))
    
}

squareWeight <- function(at, th, weight1) sum((at[at > th] - th * weight1)^2)
 

