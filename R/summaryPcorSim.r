summary.pcorSim <- function(object,...){
  # stop if the object is not a dblm object.
  if (!inherits(object, "pcorSim"))
    stop("use only with \"pcorSim\" object")
  p  <- dim(object$path)[1]
  ed <- sum(object$path)/2
  sp <- round(1 - ed/((p*p-1)/2),5)
  
  cat(gettextf("\npattern: \"%s\", Number of nodes = %s, \t Number of edges = %s, \t Sparsity = %s", object$pattern, p,ed,sp), "\n\n")  
}
