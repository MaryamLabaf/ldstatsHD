#############  equal correlation test ##########################################
equalCorrelationMatrices.test <- function(D1, D2, testStatistic = c("AS", "max", "exc"), 
                     					 testNullDist = c("asyIndep","asyDep", "np"), nite= 500, 
                     					 paired = FALSE, threshold = 2.3, exact=FALSE, excAdj = TRUE,
                     					 conf.level = 0.95, norm.approx = FALSE, ...)
{
  ## Data scaling
  D1    <- scale(D1)
  D2    <- scale(D2)
  N     <- dim(D1)[1]
  P     <- dim(D1)[2]

  if(!paired) thetaKnown <- 0
  else thetaKnown <- NULL
  
  ## Fisher Transform correlation differences 
  Test  <- ztransfCorrDiff(D1, D2, thetaKnown = thetaKnown )
  
  theta <- Test$theta
  Test  <- Test$Ts
  P2    <- length(Test)

  ## Permutations
  if(any(testNullDist=="asyDep") | any(testNullDist=="np"))
  {
    pb       <- txtProgressBar(min = 0, max = nite, style = 3)

    ## permutations
    TestsP <- array(0,dim=c(nite,P2))
    ThetaP <- array(0,dim=c(nite,P2))
    for(k in 1:nite){
       setTxtProgressBar(pb, k)
       if(exact){
        ID  <- sample(1:N)
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
       Test2  <- ztransfCorrDiff(D1e, D2e, thetaKnown = thetaKnown )

       TestsP[k,] <- (Test2$Ts)
       ThetaP[k,] <- (Test2$theta)
    }
    close(pb)
  }
   else
    TestsP <- NULL

  ## Functions call
  obj <- list()
  
  ## Tests
   if(any(testNullDist == "asyDep"))
   {
    if(any(testStatistic == "AS"))
      obj$AS$ad  <- equalCorrelationsTestSS(Test = Test, theta = theta, TestsP = TestsP, sumSquares = TRUE, dependency = TRUE, 
      										conf.level = conf.level, norm.approx = norm.approx)
    if(any(testStatistic == "max"))
      obj$MAX$ad  <- equalCorrelationsTestMax(Test = Test, theta = theta, TestsP = TestsP, dependency = TRUE, N = N, P = P, 
                        					  nite = nite, psiAdj = FALSE, thetaKnown = thetaKnown, conf.level = conf.level)
    if(any(testStatistic == "exc"))
     obj$EXC$ad  <- equalCorrelationsTestExc(Test = Test, theta = theta, TestsP = TestsP, dependency = TRUE, EX = threshold, 
										     conf.level = conf.level, excAdj = excAdj)
   }
   
   if(any(testNullDist == "asyIndep"))
   {
    if(any(testStatistic=="AS"))
      obj$AS$ai  <- equalCorrelationsTestSS(Test = Test, theta = NULL, TestsP = NULL, sumSquares = TRUE, dependency = FALSE, 
      										conf.level = conf.level, norm.approx = norm.approx)
    if(any(testStatistic == "max"))
      obj$MAX$ai  <- equalCorrelationsTestMax(Test = Test, theta = NULL, TestsP = NULL, dependency = FALSE, N = N, P = P, 
                        					  nite = 0, psiAdj = FALSE, thetaKnown = 0, conf.level = conf.level) 
    if(any(testStatistic == "exc"))
      obj$EXC$ai  <- equalCorrelationsTestExc(Test = Test, theta = NULL, TestsP = NULL, dependency = FALSE, EX = threshold, 
      									      conf.level = conf.level, excAdj = excAdj, nite = nite)
   } 
   
   if(any(testNullDist=="np"))
   {
     resAUX  <- equalCorrelationsTestPerm(Test = Test, theta = theta, TestsP = TestsP, sumSquares = TRUE, dependency = TRUE, 
      									  conf.level = conf.level, excAdj = excAdj)
     
     if(any(testStatistic=="AS"))    obj$AS$np 	<- resAUX$AS
     if(any(testStatistic=="max"))   obj$MAX$np <- resAUX$MAX      
     if(any(testStatistic=="exc"))   obj$EXC$np <- resAUX$EXC      
   }
   
  return(obj)
}


