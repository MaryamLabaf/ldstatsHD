
print.eqCorTestByRows <- function(x, ...){
  # stop if the object is not a pcorSimJoint object.
  if (!inherits(x, "eqCorTestByRows"))
    stop("use only with \"eqCorTestByRows\" object")
  if(is.null(x$paired)){
   cat(gettextf("\n \t Test for non-zero correlation matrix rows"),"\n\n")  
  }
  else{
   if(x$paired)
     callit <- "paired"
   else
     callit <- "independent"
   cat(gettextf("\n \t Test for equality of correlation matrix rows using %s data", callit),"\n\n")  
  }
  
  
  if(any(x$testStatistic == "AS"))
    cat(gettextf("number of significant rows for Tas: %s at %s conf.level, expected %s", sum(x$ciAS[,1]>0), round(x$conf.level,3), round(dim(x$ciAS)[1]*(1-x$conf.level),2)),"\n")  
  
  if(any(x$testStatistic == "max"))
    cat(gettextf("number of significant rows for Tm: %s at %s conf.level, expected %s", sum(x$ciMax[,1]>0), round(x$conf.level,3), round(dim(x$ciMax)[1]*(1-x$conf.level),2)),"\n\n")  

}
