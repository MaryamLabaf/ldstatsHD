
print.pcorSimJoint <- function(x,...){
  # stop if the object is not a pcorSimJoint object.
  if (!inherits(x, "pcorSimJoint"))
    stop("use only with \"pcorSimJoint\" object")
  p  <- dim(x$omega1)[1]
  edC <- (sum(x$omega1!=0 & x$omega2 !=0) -p)/2
  edD <- (sum(x$omega1!=0 & x$omega2 ==0))/2 + (sum(x$omega1==0 & x$omega2 !=0))/2
  spC <- round(1 - edC/((p*p-1)/2),5)
  spD <- round(1 - edD/((p*p-1)/2),5)
  
  if(x$diagCCtype =="ind")
   cat(gettextf("\nPattern: \"%s\", \t DataDepend = \"%s\" ", x$pattern,x$dataDepend), "\n")  
  else
   cat(gettextf("\nPattern: \"%s\", \t DataDepend = \"%s\", \t DiagCCtype = \"%s\" ", x$pattern, x$dataDepend,x$diagCCtype), "\n") 
  
  cat(gettextf("Number of nodes = %s, \t Common edges = %s, \t Sparsity common network = %s", p, edC, spC ), "\n")  
  cat(gettextf("Differential edges = %s, \t Sparsity differential network = %s", edD, spD ), "\n\n")  

}
