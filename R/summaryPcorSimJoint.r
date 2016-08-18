
summary.pcorSimJoint <- function(object,...){
  # stop if the object is not a pcorSimJoint object.
  if (!inherits(object, "pcorSimJoint"))
    stop("use only with \"pcorSimJoint\" object")
  p  <- dim(object$omega1)[1]
  edC <- (sum(object$omega1!=0 & object$omega2 !=0) -p)/2
  edD <- (sum(object$omega1!=0 & object$omega2 ==0))/2 + (sum(object$omega1==0 & object$omega2 !=0))/2
  spC <- round(1 - edC/((p*p-1)/2),5)
  spD <- round(1 - edD/((p*p-1)/2),5)
  
  if(object$diagCCtype =="ind")
   cat(gettextf("\nPattern: \"%s\", \t DataDepend = \"%s\" ", object$pattern,object$dataDepend), "\n")  
  else
   cat(gettextf("\nPattern: \"%s\", \t DataDepend = \"%s\", \t DiagCCtype = \"%s\" ", object$pattern, object$dataDepend,object$diagCCtype), "\n") 
  
  cat(gettextf("Number of nodes = %s, \t Common edges = %s, \t Sparsity common network = %s", p, edC, spC ), "\n")  
  cat(gettextf("Differential edges = %s, \t Sparsity differential network = %s", edD, spD ), "\n\n")  

}
