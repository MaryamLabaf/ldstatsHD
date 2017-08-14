
print.eqCorrMatTest <- function(x, ...){
  # stop if the object is not a pcorSimJoint object.
  if (!inherits(x, "eqCorrMatTest"))
    stop("use only with \"eqCorrMatTest\" object")
  if(attr(x, "identity")){
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
  for(i in 1:length(x$EXC$threshold)){
   if(any(x$testNullDist == "asyIndep"))
    cat(gettextf("asyIndep \t thr = %s, \t Texc = %s, \t pval =  %s, \t %s percent CI: %s %s ",  
    	round(x$threshold[i],dec),round(x$EXC$ai$testExc[i],dec), round(x$EXC$ai$pvalue[i],dec), round(x$conf.level * 100), round(x$EXC$ai$ci[i,1],dec) , round(x$EXC$ai$ci[i,2],dec)  ), "\n")
   if(any(x$testNullDist == "asyDep"))
    cat(gettextf("asyDep \t\t thr = %s,  \t Texc = %s, \t pval =  %s, \t %s percent CI: %s %s ",  
    	round(x$threshold[i],dec), round(x$EXC$ad$testExc[i],dec), round(x$EXC$ad$pvalue[i],dec), round(x$conf.level * 100), round(x$EXC$ad$ci[i,1],dec) , round(x$EXC$ad$ci[i,2],dec)  ), "\n")
   if(any(x$testNullDist == "np"))
    cat(gettextf("np \t\t\t thr = %s, \t Texc = %s, \t pval =  %s, \t %s percent CI: %s %s ",  
    	round(x$threshold[i],dec), round(x$EXC$np$testExc[i],dec), round(x$EXC$np$pvalue[i],dec), round(x$conf.level * 100), round(x$EXC$np$ci[i,1],dec) , round(x$EXC$np$ci[i,2],dec)  ), "\n")
  }
  cat("\n")
 }
}
