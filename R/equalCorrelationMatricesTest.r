#############  equal correlation test ##########################################
eqCorrMatTest <- function(D1, D2 = NULL, testStatistic = c("AS", "max", "exc"), 
                     	  testNullDist = c("asyIndep","asyDep", "np"), nite= 500, 
                     	  paired = FALSE, threshold = 2.3, excAdj = TRUE, exact = FALSE,
                     	  conf.level = 0.95, saddlePoint = FALSE, MINint =2, MAXint=100, ...)
{
 
  ## Checks 
  testStatistic <- controlseqCorTestByRows(D1 = D1, D2 = D2, testStatistic = testStatistic,
  										   nite = nite, paired = paired, conf.level = conf.level)
  										    
  ptestNullDist <- c("asyIndep","asyDep", "np")
  testNullDist2 <- match.arg(testNullDist, ptestNullDist, several.ok = TRUE)
  if(length(testNullDist2) <  length(testNullDist) ){
   if(length(testNullDist2) >= 1)
    warning("not all elements in ptestNullDist are well defined")
  }
  testNullDist <- testNullDist2
  
  ## functions
  if(!is.null(D2))
  {
   identity <- FALSE
   obj <- equalCorrelationMatrices.test(D1 = D1, D2 = D2, testStatistic = testStatistic,
                     					 testNullDist = testNullDist, nite = nite, 
                     					 paired = paired, threshold = threshold, excAdj = excAdj,
                     					 exact = exact, conf.level = conf.level, saddlePoint = saddlePoint,
                     					 MINint = MINint, MAXint = MAXint)
   obj$paired <- paired
  }	
  if(is.null(D2))
  {
   identity <- TRUE
   obj <- identCorrelationMatrices.test(D1 = D1, testStatistic = testStatistic,
                     					 testNullDist = testNullDist, nite = nite, 
                     					 paired = FALSE, threshold = threshold, excAdj = excAdj,
                     					 exact = exact, conf.level = conf.level)	
  }
  
 
  attr(obj, "identity")		<- identity
  obj$testStatistic <- testStatistic
  obj$testNullDist  <- testNullDist
  obj$threshold     <- threshold
  obj$conf.level    <- conf.level


  
  
  class(obj) <- "eqCorrMatTest"
   
  return(obj)
}


