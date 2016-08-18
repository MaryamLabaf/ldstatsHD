  
print.pcorSim <- function(x,...){
  # stop if the object is not a dblm object.
  if (!inherits(x, "pcorSim"))
    stop("use only with \"pcorSim\" object")
  p  <- dim(x$path)[1]
  ed <- sum(x$path)/2
  sp <- round(1 - ed/((p*p-1)/2),5)
  
  cat(gettextf("\npattern: \"%s\", Number of nodes = %s, \t Number of edges = %s, \t Sparsity = %s", x$pattern, p,ed,sp), "\n\n")  
}
