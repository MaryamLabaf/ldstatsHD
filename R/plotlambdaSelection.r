plot.lambdaSelection <- function(x, ...){
  	if (!inherits(x, "lambdaSelection")) stop("use only with \"lambdaSelection\" object")
  	
  	if(length(x$criterion)==2)
  	{
  	  ask = prod(par("mfcol")) < length(x$criterion) 
      one.fig <- prod(par("mfcol")) == 1
      if (ask) {
        oask <- devAskNewPage(TRUE)  # ask for new page
        on.exit(devAskNewPage(oask)) # exit the format plot (if true) 
      }
  	  for (i in 1:2) plot(x$lambda, x$crit.coef[[i]], ...)
  	}
  	else{
  	  	plot(x$lambda, x$crit.coef,...)
  	}  
  	 	
  			
}