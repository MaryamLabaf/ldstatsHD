
print.wfgl <- function(x, ...){
  # stop if the object is not a pcorSimJoint object.
  if (!inherits(x, "wfgl"))
    stop("use only with \"wfgl\" object")


   if(x$paired)
     callit <- "paired"
   else
     callit <- "independent"
  
   if(length(x$lambda1) == 1)
   {	   
       P	<- dim(x$path[[1]])[1]
	   edC 	<- sum(x$path[[1]] & x$path[[2]])
	   spC 	<- 1- edC/(P*(P-1))
   
	   edD1 	<- sum(x$path[[1]] & !x$path[[2]]) 
	   edD2 	<- sum(!x$path[[1]] & x$path[[2]]) 
	   edD		<- edD1 + edD2
	   spD		<- 1- edD/(P*(P-1))	   
   
	   cat(gettextf("\n \t joint partial correlation estimator using %s data", callit),"\n\n")  
	   cat(gettextf("Number of nodes = %s, \t Total number of possible edges = %s", P, P*(P-1)/2 ), "\n\n")  
	   cat(gettextf("Estimated common edges = %s, \t Sparsity estimated common network = %s", edC/2, round(spC,5) ), "\n\n")  
	   cat(gettextf("Estimated differential edges = %s, \t Sparsity estimated differential network = %s", edD/2, round(spD,5) ), "\n\n")  
	   cat(gettextf("Estimated edges for only pop.1 = %s, \t Estimated edges for only pop.2 = %s", edD1/2, edD2/2 ), "\n\n")  
	   if(x$paired&x$automLambdas) cat(gettextf("alpha2 = %s", round(x$alpha2,5)), "\n\n")  
  }
  else{	   
       P	<- dim(x$path[[1]][[1]])[1]
	   edC 	<- unlist(lapply(x$path, function(y) sum(y[[1]] & y[[2]])))
	   spC 	<- 1- edC/(P*(P-1))
      
	   edD1 	<- unlist(lapply(x$path, function(y) sum(y[[1]] & !y[[2]])))
	   edD2 	<- unlist(lapply(x$path, function(y) sum(!y[[1]] & y[[2]])))
	   edD		<- edD1 + edD2
	   spD		<- 1- edD/(P*(P-1))	   

	   cat(gettextf("\n \t joint partial correlation estimator using %s data", callit),"\n\n")  
	   cat(gettextf("lambda1 sequence of length %s", length(x$lambda1) ), "\n\n")  
	   cat(gettextf("Number of nodes = %s, \t Total number of possible edges = %s", P, P*(P-1)/2 ), "\n\n")  
	   cat(gettextf("Estimated common edges : %s -> %s, \t Sparsity estimated common network : %s -> %s",
	    min(edC/2), max(edC/2), round(min(spC),5),round(max(spC),5) ), "\n\n")  
	   cat(gettextf("Estimated differential edges : %s -> %s, \t Sparsity estimated differential network : %s -> %s", 
	   min(edD/2), max(edD/2), round(min(spD),5),round(max(spD),5) ), "\n\n") 
	   cat(gettextf("Estimated edges for only pop.1 : %s -> %s, \t Estimated edges for only pop.2 : %s -> %s", 
		 min(edD1/2), max(edD1/2), round(min(edD2),5),round(max(edD2),5) ), "\n\n") 

  }
 
}
