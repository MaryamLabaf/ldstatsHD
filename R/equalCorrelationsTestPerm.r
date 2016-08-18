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
      maxsR2   	<- apply(abs(TestsP), 1, function(x) squareWeight(x, DA, weight1))
	  testEXC 	<- squareWeight(Test, DA, weight1)
  	  Res1   	<- mean(testEXC < maxsR2)
  	  Res2		<- quantile( testEXC - maxsR2, c(1-conf.level, 1-0.0001)) 
	  c(Res1, Res2)
	}))
  }
  else
  {
	  DA 		<- EX
	  maxsR2   	<- apply(abs(TestsP), 1, function(x) squareWeight(abs(x), DA, weight1))
	  testEXC 	<- squareWeight(abs(Test), DA, weight1)
	  pvalEXC   <- mean(testEXC < maxsR2)
	  ciExc		<- quantile(testEXC - maxsR2, c(1-conf.level, 1-0.0001)) 
  }
  
  obj     <- list()
  obj$AS  <- list(testAS = mean(TestSS - SSs), pvalue = pvalSS, ci = ciSS)
  obj$MAX <- list(testMax = mean(TestMax - maxs), pvalue = pvalMax, ci = ciMax)
  obj$EXC <- list(testExc = mean(testEXC - maxsR2), pvalue = pvalEXC, ci = ciExc)
  return(obj)
}













#-log(1-pgpd2(abs(a1)[abs(a1)>EX]-EX,P1$mle[2],0,P1$mle[1]))
#(Test^2/2 +log(abs(Test)))[abs(Test)>DA] -DA

#k1 <- -log(1-pgpd2(abs(a1)[abs(a1)>EX]-EX,P1$mle[2],0,P1$mle[1]))
#k2 <- -log(2-2*pnorm((abs(a1)[abs(a1)>EX])))
#plot(abs(a1)[abs(a1)>EX],k1/k2)

#P1        <- gpd.fit(abs(rnorm(100000000)),threshold=EX)
#r1<-rnorm(1000000)
#a1 <- 1-(pgpd2(abs(r1)[abs(r1)>EX]-EX,P1$mle[2],0,P1$mle[1]))
#a2 <-  (2-2*pnorm(abs(r1)[abs(r1)>EX])) /(2-2*pnorm(EX))
#a3<-(r1^2+log(abs(r1)))[abs(r1)>EX]
#a2 <-  (2-2*pnorm(abs(r1)[abs(r1)>EX]))
# plot(abs(r1)[abs(r1)>EX],-log(a2)/a3)
 
#k1 <- -log(1-)
#k2 <- -log(1-1*pnorm((abs(r1)[abs(r1)>EX])))
#k3 <- (r1^2)[abs(r1)>EX]
#plot(abs(r1)[abs(r1)>EX],k2/k3)

#,P1$mle[2],0,P1$mle[1]))

