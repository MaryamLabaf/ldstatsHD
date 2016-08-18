identCorrelationMatrices.test  <- function(D1, testStatistic = c("AS", "max", "exc"), 
                     					 testNullDist = c("asyIndep", "asyDep", "np"), nite= 500, 
                     					 threshold = 2.3, excAdj = TRUE, conf.level = 0.95, 
                     					 norm.approx = FALSE, ...)
{
  ## Data scaling
  D1    <- scale(D1)
  N     <- dim(D1)[1]
  P     <- dim(D1)[2]

  ## Fisher Transform correlation differences 
  Test  <- ztransf(lowerTri(cor(D1))) * sqrt(N-3)
  P2    <- length(Test)

  ## Permutations
  if(any(testNullDist=="asyDep") | any(testNullDist=="np"))
  {
    pb       <- txtProgressBar(min = 0, max = nite, style = 3)

    ## Monte Carlo
    TestsP <- array(0,dim=c(nite,P2))
    for(k in 1:nite){
       setTxtProgressBar(pb, k)
       D2 		  <- mvrnorm(N, rep(0,P), diag(rep(1,P)))
       TestsP[k,] <- ztransf(lowerTri(cor(D2))) * sqrt(N-3)
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
      obj$AS$ad  <- equalCorrelationsTestSS(Test = Test, theta = NULL, TestsP = TestsP, sumSquares = TRUE, dependency = TRUE, 
      										conf.level = conf.level, norm.approx = norm.approx)
    if(any(testStatistic == "max"))
      obj$MAX$ad  <- equalCorrelationsTestMax(Test = Test, theta = NULL, TestsP = TestsP, dependency = TRUE, N = N, P = P, 
                        					  nite = nite, psiAdj = FALSE, thetaKnown = NULL, conf.level = conf.level)
    if(any(testStatistic == "exc"))
     obj$EXC$ad  <- equalCorrelationsTestExc(Test = Test, theta = NULL, TestsP = TestsP, dependency = TRUE, EX = threshold, 
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
     resAUX  <- equalCorrelationsTestPerm(Test = Test, theta = NULL, TestsP = TestsP, sumSquares = TRUE, dependency = TRUE, 
      									  conf.level = conf.level, excAdj = excAdj)
     
     if(any(testStatistic=="AS"))    obj$AS$np 	<- resAUX$AS
     if(any(testStatistic=="max"))   obj$MAX$np <- resAUX$MAX      
     if(any(testStatistic=="exc"))   obj$EXC$np <- resAUX$EXC      
   }
   
  return(obj)
}


