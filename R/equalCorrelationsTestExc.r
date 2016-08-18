#####################
equalCorrelationsTestExc <- function(Test, theta = NULL, TestsP = NULL, EX = 2.3, excAdj = TRUE,
									 dependency = TRUE, conf.level = 0.95, nite = 200, ...)
{
      pgpd2 	<- fExtremes::pgpd
      aT    	<- Test
      P2    	<- length(Test)
      weight1 	<- ifelse(excAdj, 1, 0)
      
      if(dependency)
      {

         if(length(EX)>1){
            PVALandCI  <- t(apply(as.matrix(EX),1,function(DA){
	             maxsR2   	<- apply(abs(TestsP), 1, function(x) squareWeight(x, DA, weight1))
			     maxsT 		<- squareWeight(abs(aT), DA, weight1)
            	 para 		<- fitdistr(maxsR2 + 0.000001, "gamma") 
				 Res2  		<- 1- pgamma(maxsT, para[1]$estimate[1], para[1]$estimate[2])
    	         Res3  		<- maxsT - qgamma(c(conf.level, 0.0001), para[1]$estimate[1], para[1]$estimate[2])
    	         c(Res2, Res3, maxsT - para[1]$estimate[1]/para[1]$estimate[2])
            }))
            pvalue <- PVALandCI[,1]
            ci     <- PVALandCI[,2:3]
            maxsT  <- PVALandCI[,4]
         }

         else{
             DA 		<- EX
   	     	 maxsT 		<- squareWeight(abs(aT), DA, weight1)
             maxsR2   	<- apply(abs(TestsP), 1, function(x) squareWeight(abs(x), DA, weight1))
			 para 		<- fitdistr(maxsR2 + 0.000001, "gamma") 
			 pvalue  	<- 1 - pgamma(maxsT, para[1]$estimate[1], para[1]$estimate[2])
    	     ci  		<- maxsT - qgamma(c(conf.level, 0.0001), para[1]$estimate[1], para[1]$estimate[2])
         	 maxsT  	<- maxsT - para[1]$estimate[1]/para[1]$estimate[2]
         }

      }
      else{
        if(length(EX)>1){
            	PVALandCI  <- t(apply(as.matrix(EX), 1, function(DA){
				mcarlo 		<- apply(as.matrix(1:nite),1,function(j) squareWeight(abs(rnorm(P2,0,1)), DA, weight1))
				mu0    		<- mean(mcarlo)
            	sigma2 		<- var(mcarlo)
				k	    	<- mu0^2/sigma2
				b      		<- sigma2/mu0
				Res2  		<- 1 - (pgamma(maxsT, k, 1/b))
				Res3  		<- maxsT - qgamma(c(conf.level, 0.0001), k, 1/b)
    	        c(Res2, Res3, mu0)
            }))
            pvalue <- PVALandCI[,1]
            ci     <- PVALandCI[,2:3]
            maxsT  <- maxsT - PVALandCI[,4]
         }
         else{
            DA 			<- EX
            maxsT       <- squareWeight(abs(aT), DA, weight1)
			mcarlo 		<- apply(as.matrix(1:nite),1,function(j) squareWeight(abs(rnorm(P2,0,1)), DA, weight1))
			mu0    		<- mean(mcarlo)
            sigma2 		<- var(mcarlo)
            k	    	<- mu0^2/sigma2
            b      		<- sigma2/mu0
			pvalue  	<- 1 - (pgamma(maxsT, k, 1/b))
    	    ci  		<- maxsT - qgamma(c(conf.level, 0.0001), k, 1/b)
    	    maxsT  		<- maxsT - mu0
		}
      }
      return(list(testExc = maxsT, pvalue = pvalue, ci = ci))
    
}

squareWeight <- function(at, th, weight1) sum((at[at > th] - th * weight1)^2)
 

