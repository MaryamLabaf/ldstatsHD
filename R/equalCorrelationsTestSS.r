###### Equality of correlation test: sum of squares ############################
equalCorrelationsTestSS  <- function(Test, theta = NULL, TestsP = NULL, sumSquares = TRUE, 
									 dependency = TRUE, conf.level = 0.95, norm.approx = FALSE)
{

      ## Test statistic
      P2  <- length(Test)
      if(sumSquares) TestSS <- mean(Test^2)
      else           TestSS <- var(Test)

      ## null distribution
      if(dependency){
            #NI 		<- 10000
			#means  	<- meanCors(t(TestsP)[sample(P2,min(P2,NI)),], TRUE, FALSE)
            #gamma2 	<-  2*(means[2] - NI/NI^2) * (NI/(NI-1))
            #gamma1 	<-  (means[1] - NI/NI^2) * (NI/(NI-1))

            ## moments and variance
            if(sumSquares) varBias <-  mean(apply(TestsP,2,function(x) mean(x^2)))
            else           varBias <-  mean(apply(TestsP,2,var))

            mu2 <- varBias
            mu4 <- mean(apply(TestsP^4, 1, mean))

            #if(sumSquares) VarEst  <-  max( (mu4-mu2)/P2 + (P2-1)/(P2) * gamma2,2/P2)
            #else           VarEst  <- max(mu4*(P2-2)/(P2-1)^2 - 1/(P2-1)*mu2 +(P2-2)/(P2-1) * gamma2  
            #							+ gamma1^2 + 2 * gamma1 - 2*(P2-2)/(P2-1) * gamma1,2/P2)

            if(sumSquares) VarEst <-  var(apply(TestsP,1,function(x) mean(x^2)))
            else           VarEst <-  var(apply(TestsP,1,var))

            ## pvalues approximation
            pvalueNorm      <- 1-pnorm((TestSS-varBias)/sqrt(var(apply(TestsP,1,var))),0,1)
            ciNorm			<- TestSS - qnorm(c(conf.level, 0.0001), varBias, sqrt(var(apply(TestsP,1,var))))
            #alpha           <- (varBias^2 + 2*VarEst)/VarEst
            #betaa           <- varBias*(alpha-1)
            #pvalueiGamma    <- 1-pigamma(TestSS, alpha, betaa)
			k	    		<- varBias^2/VarEst
            b       		<- VarEst/varBias
            pvalueGamma    	<- 1-pgamma(TestSS, k, 1/b)
            ciGamma			<- TestSS - qgamma(c(conf.level, 0.0001), k, 1/b)
 		    TestSS 			<- TestSS - k * b
      }
      else{
            pvalueNorm      <- 1-pnorm((TestSS-1)/sqrt(2/P2))
            ciNorm			<- TestSS - qnorm(c(conf.level, 0.0001),1, sqrt(2/P2))
            #alpha           <- (1^2 + 2*2/P2)/(2/P2)
            #betaa           <- 1*(alpha-1)
            #pvalueiGamma    <- 1-pigamma(TestSS, alpha, betaa)
			k	    		<- (P2/2)
            b       		<- 2/P2
            pvalueGamma    	<- 1-pgamma(TestSS, k, 1/b)
            ciGamma			<- TestSS - qgamma(c(conf.level, 0.0001), k, 1/b)
            TestSS 			<- TestSS - k*b
      }
      pvalue <- ifelse(norm.approx, pvalueNorm, pvalueGamma)
      if(norm.approx)       ci     <-  ciNorm
      else					ci	   <-  ciGamma
      
      return(list(testAS = TestSS, pvalue = pvalue, ci = ci))
}

