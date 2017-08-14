
###### Equality of correlation test: absolute maximum ############################
equalCorrelationsTestMax <- function(Test, theta = NULL, TestsP = NULL, dependency = TRUE, 
                                     N, P, nite, psiAdj = FALSE, thresPsi = 1.05, thetaKnown = NULL,
                                     conf.level = 0.95)
{
      pgev2 <- evd::pgev
      qgev2 <- evd::qgev
      
      ## Test statistic
      P2    <- length(Test)
      Tm    <- max(abs(Test)) 
      psi   <- NULL
      if(dependency){
         maxs  <- apply(abs(TestsP),1,max)
         GEV1  <- fgev(maxs)
         	
		 if(psiAdj){
         	D1     <- mvrnorm(N,rep(0,P),diag(rep(1,P)))
         	D2     <- mvrnorm(N,rep(0,P),diag(rep(1,P)))
         	P2     <- P*(P-1)/2
         	pb     <- txtProgressBar(min = 0, max = min(nite,P2), style = 3)
            maxsR <- apply(as.matrix(1:nite),1,function(k){
            			setTxtProgressBar(pb, k)
		            	ID  <- rbinom(N,1,0.5)
        			    D1e <- rbind(D1[ID==1,],D2[ID==0,])
			            D2e <- rbind(D2[ID==1,],D1[ID==0,])
            			a1  <-  ztransfCorrDiff(D1e, D2e, thetaKnown=thetaKnown)
			            max(abs(a1$Ts))
         	})
		 	close(pb)
         	GEV2  <- fgev(maxsR)
         
         	xi    <- (GEV2$estimate[2]-GEV1$estimate[2])/(GEV2$estimate[1]-GEV1$estimate[1])
         	psi   <- (GEV2$estimate[2]/GEV1$estimate[2])^(-1/xi)
         	if(psi > thresPsi)
         	{
              pnM   <- pnorm(Tm,  0,1)
              pnM2  <- pnorm(-Tm, 0,1)
              pval  <- (pnM - pnM2)^P2
              psi   <- 1
              ci    <- Tm - (sqrt(2 * log(2 * P2)) -(log(log(2*P2)) + log(pi*4*log(2,2)))/(2 * sqrt(2 * log(2 * P2))) -
          				log(- log(c(conf.level, 0.0001)))/sqrt(2 * log(2 * P2)))
          	  Tm    <-  Tm - (sqrt(2 * log(2 * P2)) -(log(log(2*P2)) + log(pi*4*log(2,2)))/(2 * sqrt(2 * log(2 * P2)))) 
         	}
        	else{
        	 pval   <- pgev2(Tm, GEV1$estimate[1], GEV1$estimate[2], (GEV1$estimate[3]+GEV2$estimate[3])/2)
        	 xi 	<- (GEV1$estimate[3]+GEV2$estimate[3])/2
        	 ci     <- Tm - qgev2(c(conf.level, 0.0001), GEV1$estimate[1], GEV1$estimate[2], xi)
	         Tm     <- Tm -  GEV1$estimate[1] + GEV1$estimate[2] * (gamma(1-xi) - 1)/xi
         	}
         }
         else{
          pval  <- pgev2(Tm, GEV1$estimate[1], GEV1$estimate[2], GEV1$estimate[3])
          ci    <- Tm - qgev2(c(conf.level, 0.0001), GEV1$estimate[1], GEV1$estimate[2], GEV1$estimate[3])
          Tm    <- Tm -  GEV1$estimate[1] + GEV1$estimate[2] * (gamma(1-GEV1$estimate[3]) - 1)/GEV1$estimate[3]
		 }
     }
     else{
      pnM   <- pnorm(Tm,0,1)
      pnM2  <- pnorm(-Tm,0,1)
      pval  <- (pnM - pnM2)^P2
      psi 	<- 1
      ci    <- Tm - (sqrt(2 * log(2 * P2)) -(log(log(2*P2)) + log(pi*4*log(2,2)))/(2 * sqrt(2 * log(2 * P2))) -
          				log(- log(c(conf.level, 0.0001)))/sqrt(2 * log(2 * P2)))
      Tm    <-  Tm - (sqrt(2 * log(2 * P2)) -(log(log(2*P2)) + log(pi*4*log(2,2)))/(2 * sqrt(2 * log(2 * P2)))) 
          				
     }
     
     return(list(testMax = Tm, pvalue = 1 - pval, psi = psi, ci = ci))
}
