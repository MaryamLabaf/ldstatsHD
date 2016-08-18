
print.lambdaSelection <- function(x, ...){
  # stop if the object is not a pcorSimJoint object.
  if (!inherits(x, "lambdaSelection"))
    stop("use only with \"lambdaSelection\" object")

   
   if(length(x$criterion)==2)
       cat(gettextf("\n \t lambda selection by optimizing %s / %s risk function", x$criterion[1], x$criterion[2]),"\n\n")  
    else
    {
     if(x$criterion == "A-MSE")
     { 
      callit <- attr(x, "generator")
      cat(gettextf("\n \t lambda selection by optimizing %s risk function with %s generator ", x$criterion, callit),"\n\n")  
     }
     else
       cat(gettextf("\n \t lambda selection by optimizing %s risk function", x$criterion),"\n\n")  
    }
    
   P  <- dim(attr(x, "bestpath"))[1] 
   if(length(x$criterion)==2)
   {
    SP1 <- sum(attr(x, "bestpath"))/(P*(P-1))
    SP2 <- sum(attr(x, "bestpath2"))/(P*(P-1))
    cat(gettextf("optimal lambda  %s = %s, \t Sparsity graph structure = %s", x$criterion[1], round(x$opt.lambda[1],4), round(1-SP1,4)), "\n\n")
    cat(gettextf("optimal lambda  %s = %s, \t Sparsity graph structure = %s", x$criterion[2], round(x$opt.lambda[2],4), round(1-SP2,4)), "\n\n")
   }
   else{ 
    SP <- sum(attr(x, "bestpath"))/(P*(P-1))
    cat(gettextf("optimal lambda = %s, \t Sparsity graph structure = %s", round(x$opt.lambda[1],4), round(1-SP,4)), "\n\n")
   }  
   
 
}
