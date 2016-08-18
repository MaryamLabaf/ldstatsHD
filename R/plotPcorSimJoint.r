plot.pcorSimJoint <- function(x, minn = 0, col = c("blue","red","green"), vertex.size=3,  edgesThickness = FALSE, ...){
  if (!inherits(x, "pcorSimJoint")) stop("use only with \"pcorSimJoint\" object")
  
  if(edgesThickness)
  {
   C1 <- abs(x$omega1)
   C2 <- abs(x$omega2)
   pcor2jointGraph(x$omega1, x$omega2, minn = minn, col = col, vertex.size = vertex.size, C1 = C1, C2=C2, edgesThicknees=TRUE, ...)
  }
  else
   pcor2jointGraph(x$omega1, x$omega2, minn = minn, col = col, vertex.size = vertex.size, edgesThicknees = FALSE, ...)
     
}
