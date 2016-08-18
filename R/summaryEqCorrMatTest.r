
summary.eqCorrMatTest <- function(object, ...){
  x <- object
  # stop if the object is not a pcorSimJoint object.
  if (!inherits(x, "eqCorrMatTest"))
    stop("use only with \"eqCorrMatTest\" object")
  if(is.null(x$paired)){
   cat(gettextf("\n \t Test for non-Identity correlation matrix"),"\n\n")  
  }
  else{
   if(x$paired)
     callit <- "paired"
   else
     callit <- "independent"
   cat(gettextf("\n \t Test for equality of two correlation matrices using %s data", callit),"\n\n")  
  }
  
  dec <- 3
  
  if(any(x$testStatistic == "AS")){
   if(any(x$testNullDist == "asyIndep"))
    cat(gettextf("asyIndep \t Tas = %s, \t pval =  %s, \t %s percent CI: %s %s ",  
    	round(x$AS$ai$testAS,dec), round(x$AS$ai$pvalue,dec), round(x$conf.level * 100), round(x$AS$ai$ci[1],dec) , round(x$AS$ai$ci[2],dec) ),  "\n")
   if(any(x$testNullDist == "asyDep"))
    cat(gettextf("asyDep \t\t Tas = %s, \t pval =  %s, \t %s percent CI: %s %s ",  
    	round(x$AS$ad$testAS,dec), round(x$AS$ad$pvalue,dec), round(x$conf.level * 100), round(x$AS$ad$ci[1],dec) , round(x$AS$ad$ci[2],dec) ),  "\n")
   if(any(x$testNullDist == "np"))
    cat(gettextf("np \t\t\t Tas = %s, \t pval =  %s, \t %s percent CI: %s %s ",  
    	round(x$AS$np$testAS,dec), round(x$AS$np$pvalue,dec), round(x$conf.level * 100), round(x$AS$np$ci[1],dec) , round(x$AS$np$ci[2],dec) ), "\n")
   cat("\n")
  }
  if(any(x$testStatistic == "max")){
   if(any(x$testNullDist == "asyIndep"))
    cat(gettextf("asyIndep \t Tm = %s, \t pval =  %s, \t %s percent CI: %s %s ",  
    	round(x$MAX$ai$testMax,dec), round(x$MAX$ai$pvalue,dec), round(x$conf.level * 100), round(x$MAX$ai$ci[1],dec) , round(x$MAX$ai$ci[2],dec) ), "\n")
   if(any(x$testNullDist == "asyDep"))
    cat(gettextf("asyDep \t\t Tm = %s, \t pval =  %s, \t %s percent CI: %s %s ",  
    	round(x$MAX$ad$testMax,dec), round(x$MAX$ad$pvalue,dec), round(x$conf.level * 100), round(x$MAX$ad$ci[1],dec) , round(x$MAX$ad$ci[2],dec) ), "\n")
   if(any(x$testNullDist == "np"))
    cat(gettextf("np \t\t\t Tm = %s, \t pval =  %s, \t %s percent CI: %s %s ",  
    	round(x$MAX$np$testMax,dec), round(x$MAX$np$pvalue,dec), round(x$conf.level * 100), round(x$MAX$np$ci[1],dec) , round(x$MAX$np$ci[2],dec) ), "\n")   
  cat("\n")
 }
 if(any(x$testStatistic == "exc")){
   if(any(x$testNullDist == "asyIndep"))
    cat(gettextf("asyIndep \t Texc = %s, \t pval =  %s, \t %s percent CI: %s %s ",  
    	round(x$EXC$ai$testExc,dec), round(x$EXC$ai$pvalue,dec), round(x$conf.level * 100), round(x$EXC$ai$ci[1],dec) , round(x$EXC$ad$ci[2],dec)  ), "\n")
   if(any(x$testNullDist == "asyDep"))
    cat(gettextf("asyDep \t\t Texc = %s, \t pval =  %s, \t %s percent CI: %s %s ",  
    	round(x$EXC$ad$testExc,dec), round(x$EXC$ad$pvalue,dec), round(x$conf.level * 100), round(x$EXC$ad$ci[1],dec) , round(x$EXC$ad$ci[2],dec)  ), "\n")
   if(any(x$testNullDist == "np"))
    cat(gettextf("np \t\t\t Texc = %s, \t pval =  %s, \t %s percent CI: %s %s ",  
    	round(x$EXC$np$testExc,dec), round(x$EXC$np$pvalue,dec), round(x$conf.level * 100), round(x$EXC$np$ci[1],dec) , round(x$EXC$np$ci[2],dec)  ), "\n")
  cat("\n")
 }
 
}
