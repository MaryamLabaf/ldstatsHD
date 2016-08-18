
print.eqCorrMatTest <- function(x, ...){
  # stop if the object is not a pcorSimJoint object.
  if (!inherits(x, "eqCorrMatTest"))
    stop("use only with \"eqCorrMatTest\" object")
  if(attr(x, "identity")){
   cat(gettextf("\n Test for non-Identity correlation matrix"),"\n\n")  
  }
  else{
   if(x$paired)
     callit <- "paired"
   else
     callit <- "independent"
   cint  <- round(x$conf.level * 100)
   cat(gettextf("\n Test for equality of two correlation matrices using %s data, %s percent confidence interval (ci): ", callit, cint),"\n\n")  
  }
  
  dec <- 3
  c1 <- numeric()
  teststat <- character()
  testnull <- character()
  
  if(any(x$testStatistic == "AS")){
   if(any(x$testNullDist == "asyIndep"))
   {
    	c1 <- rbind(c1,as.numeric(c(round(x$AS$ai$testAS,dec), round(x$AS$ai$pvalue,dec),  round(x$AS$ai$ci[1],dec) , round(x$AS$ai$ci[2],dec) )))
        teststat <- c(teststat,"AS")
        testnull <- c(testnull,"asyIndep")
   }
   if(any(x$testNullDist == "asyDep"))
   {
    	c1 <- rbind(c1,as.numeric(c(round(x$AS$ad$testAS,dec), round(x$AS$ad$pvalue,dec), round(x$AS$ad$ci[1],dec) , round(x$AS$ad$ci[2],dec) )))
        teststat <- c(teststat,"AS")
        testnull <- c(testnull,"asyDep")
   }
   if(any(x$testNullDist == "np"))
   {
    	c1 <- rbind(c1,as.numeric(c(round(x$AS$np$testAS,dec), round(x$AS$np$pvalue,dec), round(x$AS$np$ci[1],dec) , round(x$AS$np$ci[2],dec) )))
        teststat <- c(teststat,"AS")
        testnull <- c(testnull,"np")
   }
  }
  if(any(x$testStatistic == "max")){
   if(any(x$testNullDist == "asyIndep"))
   {
    	c1 <- rbind(c1,as.numeric(c(round(x$MAX$ai$testMax,dec), round(x$MAX$ai$pvalue,dec), round(x$MAX$ai$ci[1],dec) , round(x$MAX$ai$ci[2],dec) )))
        teststat <- c(teststat,"MAX")
        testnull <- c(testnull,"asyIndep")
   }
   if(any(x$testNullDist == "asyDep"))
   {
    	c1 <- rbind(c1,as.numeric(c(round(x$MAX$ad$testMax,dec), round(x$MAX$ad$pvalue,dec),round(x$MAX$ad$ci[1],dec) , round(x$MAX$ad$ci[2],dec) )))
        teststat <- c(teststat,"MAX")
        testnull <- c(testnull,"asyDep")
   }
   if(any(x$testNullDist == "np"))
   {
    	c1 <- rbind(c1,as.numeric(c(round(x$MAX$np$testMax,dec), round(x$MAX$np$pvalue,dec), round(x$MAX$np$ci[1],dec) , round(x$MAX$np$ci[2],dec) )))
        teststat <- c(teststat,"MAX")
        testnull <- c(testnull,"np")
   }
  }
 if(any(x$testStatistic == "exc")){
   if(any(x$testNullDist == "asyIndep"))
   {
    	c1 <- rbind(c1,as.numeric(c(round(x$EXC$ai$testExc,dec), round(x$EXC$ai$pvalue,dec),round(x$EXC$ai$ci[1],dec) , round(x$EXC$ai$ci[2],dec) )))
        teststat <- c(teststat,"EXC")
        testnull <- c(testnull,"asyIndep")
   }
   if(any(x$testNullDist == "asyDep"))
   {
    	c1 <- rbind(c1,as.numeric(c(round(x$EXC$ad$testExc,dec), round(x$EXC$ad$pvalue,dec),  round(x$EXC$ad$ci[1],dec) , round(x$EXC$ad$ci[2],dec) )))
        teststat <- c(teststat,"EXC")
        testnull <- c(testnull,"asyDep")
   }
   if(any(x$testNullDist == "np"))
   {
    	c1 <- rbind(c1,as.numeric(c(round(x$EXC$np$testExc,dec), round(x$EXC$np$pvalue,dec), round(x$EXC$np$ci[1],dec) , round(x$EXC$np$ci[2],dec) )))
        teststat <- c(teststat,"EXC")
        testnull <- c(testnull,"np")
   }
 }
 
 print(data.frame(test = teststat, testnull = testnull, stat = c1[,1],  p.value = c1[,2], lower.ci = c1[,3], upper.ci = c1[,4]))
 cat("\n")
}
