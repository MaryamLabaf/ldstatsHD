
print.wfrl <- function(x, ...){
  # stop if the object is not a pcorSimJoint object.
  if (!inherits(x, "wfrl"))
    stop("use only with \"wfrl\" object")


   if(x$paired)
     callit <- "paired"
   else
     callit <- "independent"
  
   if(length(x$lambda1) == 1)
   {	   
       P	<- dim(x$path[[1]])[1]
   	   Q	<- dim(x$path[[1]])[2]
	   edC 	<- sum(x$path[[1]] & x$path[[2]])
	   spC 	<- 1- edC/(Q*P)
   
	   edD1 	<- sum(x$path[[1]] & !x$path[[2]]) 
	   edD2 	<- sum(!x$path[[1]] & x$path[[2]]) 
	   edD		<- edD1 + edD2
	   spD		<- 1- edD/(Q*P)	      
	   
	   cat(gettextf("\n \t joint regression coefficients estimator using %s data", callit),"\n\n")  
       cat(gettextf("Number of response variables = %s, \t Number of explanatory variables = %s \t Number of possible edges = %s", Q, P, Q*P), "\n\n")  
	   cat(gettextf("Estimated common edges = %s, \t Sparsity estimated common network = %s", edC, round(spC,5) ), "\n\n")  
	   cat(gettextf("Estimated differential edges = %s, \t Sparsity estimated differential network = %s", edD, round(spD,5) ), "\n\n")  
	   cat(gettextf("Estimated edges for only pop.1 = %s, \t Estimated edges for only pop.2 = %s", edD1, edD2 ), "\n\n")  
  }
  else{	   
       P	<- dim(x$path[[1]][[1]])[1]
   	   Q	<- dim(x$path[[1]][[1]])[2]

	   edC 	<- unlist(lapply(x$path, function(y) sum(y[[1]] & y[[2]])))
	   spC 	<- 1- edC/(Q*P)
      
	   edD1 	<- unlist(lapply(x$path, function(y) sum(y[[1]] & !y[[2]])))
	   edD2 	<- unlist(lapply(x$path, function(y) sum(!y[[1]] & y[[2]])))
	   edD		<- edD1 + edD2
	   spD		<- 1- edD/(Q*P)	 

	   cat(gettextf("\n \t joint regression coefficients estimator using %s data", callit),"\n\n")  
	   cat(gettextf("lambda1 sequence of length %s", length(x$lambda1) ), "\n\n")  
      cat(gettextf("Number of response variables = %s, \t Number of explanatory variables = %s \t Number of possible edges = %s", Q, P, Q*P), "\n\n")  
	   cat(gettextf("Estimated common edges : %s -> %s, \t Sparsity estimated common network : %s -> %s",
	    min(edC), max(edC), round(min(spC),5),round(max(spC),5) ), "\n\n")  
	   cat(gettextf("Estimated differential edges : %s -> %s, \t Sparsity estimated differential network : %s -> %s", 
	   min(edD), max(edD), round(min(spD),5),round(max(spD),5) ), "\n\n") 
	   cat(gettextf("Estimated edges for only pop.1 : %s -> %s, \t Estimated edges for only pop.2 : %s -> %s", 
		 min(edD1), max(edD1), round(min(edD2),5),round(max(edD2),5) ), "\n\n") 


  }
 
}


