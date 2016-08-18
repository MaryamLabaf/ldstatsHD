plot.wfgl <- function(x, minn = 0, col = c("blue","red","green"), vertex.size = 3,  edgesThickness = FALSE, zoomThick = 10, ...){
  if (!inherits(x, "wfgl")) stop("use only with \"wfgl\" object")
  
  if(length(x$lambda1)==1)
  {
	  if(edgesThickness)
	  {
	   C1 <- abs(x$omega[[1]])
	   C2 <- abs(x$omega[[2]])
	   pcor2jointGraph(x$path[[1]], x$path[[2]], minn = minn, col = col, vertex.size = vertex.size, C1 = C1, C2 = C2, 
					   edgesThicknees=TRUE, zoomThick = zoomThick, ...)
	  }
	  else
	   pcor2jointGraph(x$path[[1]], x$path[[2]], minn = minn, col = col, vertex.size = vertex.size, edgesThicknees = FALSE, ...)
  }
  else
  {
	     ask = prod(par("mfcol")) < length(x$lambda1) 
      one.fig <- prod(par("mfcol")) == 1
     if (ask) {
        oask <- devAskNewPage(TRUE)  # ask for new page
        on.exit(devAskNewPage(oask)) # exit the format plot (if true) 
     }

     for (i in 1:length(x$lambda1)){
    
	  if(edgesThickness)
	  {
	   C1 <- abs(x$omega[[i]][[1]])
	   C2 <- abs(x$omega[[i]][[2]])
	   pcor2jointGraph(x$path[[i]][[1]], x$path[[i]][[2]], minn = minn, col = col, vertex.size = vertex.size, C1 = C1, C2 = C2, 
					   edgesThicknees=TRUE, zoomThick = zoomThick, ...)
	  }
	  else
	   pcor2jointGraph(x$path[[i]][[1]], x$path[[i]][[2]], minn = minn, col = col, vertex.size = vertex.size, edgesThicknees = FALSE, ...)
      
	}  
  
  }      

}