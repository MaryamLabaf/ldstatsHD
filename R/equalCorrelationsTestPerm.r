###### Equality of correlation test: Permutation Tests ############################
equalCorrelationsTestPerm  <- function(Test, theta, TestsP, EX = 2.3, sumSquares=TRUE, excAdj = TRUE, dependency=TRUE,
										conf.level = 0.95)
{
  pgpd2 	<- fExtremes::pgpd
  P2      	<- length(Test)

  ##SS
  TestSS  	<- ifelse(sumSquares, mean(Test^2), var(Test))
  SSs     	<- apply((TestsP)^2,1,mean)
  pvalSS 	<- 1- sum(TestSS>SSs)/length(SSs)
  ciSS		<- quantile(TestSS - SSs, c(1-conf.level, 1-0.0001)) 
   
  ##MAX
  TestMax 	<- max(abs(Test))
  maxs  	<- apply(abs(TestsP),1,max)
  pvalMax 	<- 1- sum(TestMax>maxs)/length(maxs)
  ciMax		<- quantile(TestMax - maxs, c(1-conf.level, 1-0.0001)) 
 
  ##EXC
  weight1 	<- ifelse(excAdj, 1, 0)
  if(length(EX)>1){
    PVALandCI  <- t(apply(as.matrix(EX), 1, function(DA){
	  maxsR2   	<- apply(abs(TestsP), 1, function(x) squareWeight(abs(x), DA, weight1))
	  testEXC 	<- squareWeight(abs(Test), DA, weight1)
	  Res1      <- mean(testEXC < maxsR2)
  	  Res2		<- quantile( testEXC - maxsR2, c(1-conf.level, 1-0.0001)) 
	  ve		<- var(maxsR2) 
	  c(Res1, Res2, mean(testEXC - maxsR2))
	}))

	pvalEXC <- PVALandCI[,1]
    saf  <- t(apply(as.matrix(EX), 1, function(DA){
	  maxsR2   	<- apply(abs(TestsP), 1, function(x) squareWeight(abs(x), DA, weight1))
	  testEXC 	<- squareWeight(abs(Test), DA, weight1)
	  c((maxsR2-mean(maxsR2))/sd(maxsR2),(testEXC-mean(maxsR2))/sd(maxsR2))
	}))
	pvalEXC <- c(1-mean(apply(saf[,1:dim(TestsP)[1]],2,max) <= max(saf[,dim(TestsP)[1]+1])),pvalEXC)

	ciExc   <- 1#PVALandCI[,c(2,3)]
	TESTex  <- saf[,dim(TestsP)[1]+1]#PVALandCI[,4]

  }
  
  else
  {
	  DA 		<- EX
	  maxsR2   	<- apply(abs(TestsP), 1, function(x) squareWeight(abs(x), DA, weight1))
	  testEXC 	<- squareWeight(abs(Test), DA, weight1)
	  pvalEXC   <- mean(testEXC < maxsR2)
	  ciExc		<- t(as.matrix(quantile(testEXC - maxsR2, c(1-conf.level, 1-0.0001)) ))

#	  maxsR2   	<- apply(abs(TestsP), 1, function(x) squareWeight(abs(x), DA, 1-weight1))
#	  testEXC 	<- squareWeight(abs(Test), DA,  1-weight1)
#	  pval2   <- mean(testEXC < maxsR2)
	  TESTex  <-  mean(testEXC - maxsR2)

  }
  
  obj     <- list()
  obj$AS  <- list(testAS = mean(TestSS - SSs), pvalue = pvalSS, ci = ciSS)
  obj$MAX <- list(testMax = mean(TestMax - maxs), pvalue = pvalMax, ci = ciMax)
  obj$EXC <- list(testExc = TESTex, pvalue = pvalEXC, ci = ciExc)#, pvalue2 = pval2)
  return(obj)
}

