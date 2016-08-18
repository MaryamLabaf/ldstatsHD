eqCorrelationTest.byRows <- function(D1, D2 = NULL, testStatistic = c("AS", "max"), nite = 200, 
							paired = FALSE, exact = TRUE, subMatComp = FALSE, iniP = 1, 
							finP = NULL, conf.level = 0.95)
{
  P  <- dim(D1)[2]
  N  <- dim(D1)[1]
  D1 <- scale(D1)
  D2 <- scale(D2)
   
  if(!paired) thetaKnown <- 0
  else thetaKnown <- NULL
    
  ## diag crosCorr
  diagCor <- apply(D1*D2,2,mean)/(N-1)*N
  
  if(subMatComp)
  {
   ## first to last variables to test
   j  <- iniP - 1
   if(is.null(finP)) finP <- P
   pvalSS   <- numeric(finP-j)
   pvalMax  <- numeric(finP-j)
   SStestT  <- numeric(finP-j)
   MaxTestT <- numeric(finP-j)
   ciMaxT   <- array(0, dim=c(finP-j,2))
   ciSST    <- array(0, dim=c(finP-j,2))
   pb        <- txtProgressBar(min = 0, max = nite * ceiling((finP-j)/20), style = 3)
   while(j < finP){
     if(N %% 2 !=0 & exact) 
       unfort <- sort(sample(1:N, N-1))
      else 
       unfort <- 1:N
        
      steps   <- min(20,finP-j)
      row     <- c(1:steps + j)
      ulr     <- 1:length(row) 
      testa   <- ztransfCorrDiffsub(D1[unfort,], D2[unfort,], row = row, 
      					thetaKnown = thetaKnown, diagCor = diagCor)
      test2   <- testa$Ts
      
      if(any(testStatistic == "AS")) 
        SStest  <- apply(as.matrix(1:steps),1,function(i) mean(test2[(i-1)*(P-1) + 1:(P-1)]^2))
      if(any(testStatistic == "max")) 
        MaxTest <- apply(as.matrix(1:steps),1,function(i) max(abs(test2[(i-1)*(P-1) + 1:(P-1)])))
 
      pvalSSaux  <- array(0,dim=c(nite,steps))
      pvalMaxaux <-  array(0,dim=c(nite,steps))
      
      for (ite in 1:nite){
        setTxtProgressBar(pb, ite + j/20*nite)    
        if(exact & paired){
         ID <- sample(1:N)
         ID1 <- ID[1:floor(N/2)]
         ID2 <- ID[1:floor(N/2) +floor(N/2)]
        }
        else{
         ID    <- rbinom(N,1,0.5)
         ID1   <- which(ID==1)
         ID2   <- which(ID==0)
        }
       
        D1e   <- scale(rbind(D1[ID1,],D2[ID2,]))
        D2e   <- scale(rbind(D2[ID1,],D1[ID2,]))
        
        diagCor2 <- c(apply(D1e[,1:steps+j] * D2e[,1:steps+j],2,mean)/(N-1)*N, diagCor)
             
        testa2   <- ztransfCorrDiffsub( cbind(D1e[,1:steps+j], D1) , cbind(D2e[,1:steps+j], D2), row = 1:steps, 
                               thetaKnown= thetaKnown, diagCor = diagCor2)
        theta2 <-  testa2$theta
        if(!paired) theta2 <- rep(theta2,length(testa2$Ts))  
        theta2 <-  theta2[- as.numeric(apply(as.matrix(c(0:(steps-1))),1,function(j) 
            					j*(length(diagCor2)-1) + 1:(steps-1)))]
        
        corDif <- ztransf(cor(D1e[,1:steps+j], D1[unfort,])) - ztransf(cor(D2e[,1:steps+j], D2[unfort,]))
        corDif <- as.numeric(t(sqrt(N-3) * corDif))[- ((ulr-1)*P +row)]/sqrt(2-2*theta2)[- ((ulr-1)*P +row)]
      
        if(any(testStatistic=="AS")) 
         pvalSSaux[ite,]   <- apply(as.matrix(1:steps),1,function(i) mean(corDif[(i-1)*(P-1) + 1:(P-1)]^2))
        if(any(testStatistic=="max")) 
         pvalMaxaux[ite,]  <- apply(as.matrix(1:steps),1,function(i) max(abs(corDif[(i-1)*(P-1) + 1:(P-1)])))         
     }
  
     if(any(testStatistic=="AS")) 
     {
      pvalSS[1:steps+j]   <- apply(apply(pvalSSaux, 1, function(x) x - SStest >0), 1, mean, na.rm=TRUE)
      SStestT[1:steps+j]  <- SStest
      ciSST[1:steps+j,]   <- t(apply(apply(pvalSSaux, 1, function(x) SStest - x), 1, quantile, c((1-conf.level), 1),na.rm =TRUE))
     }
     if(any(testStatistic=="max")) 
     {
      pvalMax[1:steps+j]  <- apply(apply(pvalMaxaux, 1, function(x) x - MaxTest >0), 1, mean, na.rm=TRUE)
      MaxTestT[1:steps+j] <- MaxTest     
      ciMaxT[1:steps+j,]  <- t(apply(apply(pvalMaxaux, 1, function(x) MaxTest - x), 1, quantile, c((1-conf.level), 1),na.rm =TRUE))
     }  
     j <- j + steps
   }
   close(pb)
   MaxTest <- MaxTestT
   SStest  <- SStestT
  }
  else{
     pb        <- txtProgressBar(min = 0, max = nite, style = 3)
     if(N %% 2 != 0 & exact) 
       unfort <- sort(sample(1:N,N-1))
     else 
       unfort <- 1:N

     test   <- ztransfCorrDiff(D1[unfort,], D2[unfort,], thetaKnown = thetaKnown)
     Theta 	<- array(0,dim=c(P,P))
	 Theta[lower.tri(Theta)] <- test$theta
	 Theta <- Theta + t(Theta)	
     
     sumSq 	<- array(0, dim = c(nite,P) )
     maxT 	<- array(0, dim = c(nite,P) )
     for(id in 1:nite){
       setTxtProgressBar(pb, id)    
       if(exact & paired){
        ID 	<- sample(1:N)
        ID1 <- ID[1:floor(N/2)]
        ID2 <- ID[1:floor(N/2) +floor(N/2)]
       }
       else{
        ID    <- rbinom(N,1,0.5)
        ID1   <- which(ID==1)
        ID2   <- which(ID==0)
       }
       D1e    <- rbind(D1[ID1,],D2[ID2,])
       D2e    <- rbind(D2[ID1,],D1[ID2,])
       corDif <- ztransf(cor(D1e,D1[unfort,])) - ztransf(cor(D2e,D2[unfort,]))
       corDif <- sqrt(N-3) * corDif/sqrt(2-2*Theta)
       diag(corDif) <- 0 
       if(any(testStatistic == "AS")) 
        sumSq[id,]   <-   apply(corDif^2, 1, sum)/(P-1) 
       if(any(testStatistic == "max")) 
        maxT[id,]    <-   apply(abs(corDif), 1, max)
     }
       
     AUX 	<- array(0,dim=c(P,P))
	 AUX[lower.tri(AUX)] <- test$Ts
	 AUX     <- AUX + t(AUX)	
     SStest  <- apply(AUX^2,1,sum)/(P-1)
     MaxTest <-apply(abs(AUX),1,max)
     
     if(any(testStatistic == "AS")) 
     {
      pvalSS  <- apply( as.matrix(1:P), 1,function(k) mean(SStest[k] < sumSq[,k]))
      ciSST   <- t(apply(apply(sumSq, 1, function(x) SStest - x), 1, quantile, c((1-conf.level), 1)))
     }
     if(any(testStatistic == "max")) 
     {
      pvalMax <- apply( as.matrix(1:P), 1,function(k) mean(MaxTest[k] < maxT[,k]))
      ciMaxT   <- t(apply(apply(maxT, 1, function(x) MaxTest - x), 1, quantile, c((1-conf.level), 1)))
     }
     close(pb)
  }
  
   obj 			<- list()
   obj$Maxtest 	<- NULL
   obj$pvalMax 	<- NULL
   obj$ciMax   	<- NULL
   obj$AStest 	<- NULL
   obj$pvalAS 	<- NULL
   obj$ciAS   	<- NULL
  if(any(testStatistic == "AS")) 
  {
   obj$AStest <- SStest
   obj$pvalAS <- pvalSS
   obj$ciAS   <- ciSST
  }
  if(any(testStatistic == "max")) 
  {
   obj$Maxtest <- MaxTest
   obj$pvalMax <- pvalMax
   obj$ciMax   <- ciMaxT
  }
  obj$paired 		<- paired
  obj$testStatistic <- testStatistic
  obj$conf.level 	<- conf.level
  
  return(obj)
}