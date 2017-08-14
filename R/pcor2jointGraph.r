### Graphical representation of two similar networks ###
pcor2jointGraph <- function(M1, M2, minn=0, col=c("blue","red","green"), vertex.size=2, C1=NULL,C2=NULL, 
							edgesThicknees=FALSE, zoomThick = 10, ...){
  A1 <- (M1!=0)*1
  A2 <- (M2!=0)*1 
  diag(A1) <- 0
  diag(A2) <- 0
  Aj <- (M1!=0|M2!=0)*1
  diag(Aj) <- 0
  MAJ<- Matrix(Aj)
  AJ2 <- MAJ%*%MAJ
  wh1 <- which(apply(Aj,1,sum)>0&apply(AJ2,1,sum)>minn)
  Aj2 <- Aj[wh1,wh1]
  wh11 <- which(apply(Aj2,1,sum)>0)
  Aj2 <- Aj2[wh11,wh11]
  
  gr <- graph.adjacency(Aj2,mode="undirected")
  NAMES <- get.data.frame(gr, what="edg")
  INwhich <- apply(NAMES,1,function(x){
    wh2 <- wh1[as.numeric(x)]
    c(M1[wh2[1],wh2[2]]!=0,M2[wh2[1],wh2[2]]!=0)
    })
  COL   <- rep(col[1],dim(NAMES)[1])
  COL[INwhich[1,]&!INwhich[2,]] <- col[2]
  COL[!INwhich[1,]&INwhich[2,]] <- col[3]
  if(edgesThicknees)
  {
    INwhich2 <- apply(NAMES,1,function(x){
     wh2 <- wh1[as.numeric(x)]
     abs(c(C1[wh2[1],wh2[2]],C2[wh2[1],wh2[2]]))
    })
   thick <- rep(1,dim(NAMES)[1])
   thick[INwhich[1,]&!INwhich[2,]] <- abs(INwhich2[1,INwhich[1,]&!INwhich[2,]]-INwhich2[2,INwhich[1,]&!INwhich[2,]])
   thick[!INwhich[1,]&INwhich[2,]] <- abs(INwhich2[1,!INwhich[1,]&INwhich[2,]]-INwhich2[2,!INwhich[1,]&INwhich[2,]])
   thick[INwhich[1,]&INwhich[2,]]  <- INwhich2[1,INwhich[1,]&INwhich[2,]]

   plot(gr,vertex.size=vertex.size,vertex.label=NA,edge.color= COL, edge.width = thick * zoomThick, ...)
  }
  else
   plot(gr,vertex.size=vertex.size,vertex.label=NA,edge.color= COL,...)
   
}
