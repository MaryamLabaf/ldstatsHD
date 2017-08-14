pvaluesTriangleWeakestEdge <- function(A1, A2, alpha=0.01, D1, D2, equalBetweenCond = TRUE)
{
          ## Find Triangles
          inda  <- array(rep(1:3,3),dim=c(3,3))
          ind1  <- lowerTri(t(inda))
          ind2  <- lowerTri(inda)
          N <- dim(D1)[1]
          diag(A1)          <-0
          diag(A2)          <-0
          A1.new            <- A1
          A2.new            <- A2
          edgesToConsider   <- vector(2,mode="list")
          pvaluesToConsider <- vector(2,mode="list")

          TRI1r           <- findTriangles(A1)
          TRI2r           <- findTriangles(A2)
		 print(dim(TRI1r))
		 print(dim(TRI2r))
		 
          ## pvalues first matrix
          cor2EX <- FALSE
          if(class(TRI1r)=="matrix"){
             if(dim(TRI1r)[1]>0){
                cor2 <- apply(as.matrix(1:dim(TRI1r)[1]),1,function(PL){
                  print(PL)
                  cc1 <- TRI1r[PL,]
                  aa <- try(abs(cov2cor(solve(cov2cor(cor(D1[,cc1]))))))
                  if(class(aa)!="try-error"){
                   aa2 <- which.min(lowerTri(aa))
                   aa  <- min(aa)
                   return(c(c(cc1[ind1[aa2]],cc1[ind2[aa2]]),pnorm(ztransf(aa)*sqrt(N-5))))
   				  }
   				  else
   				  {
   				   aa2 <- 1
                   aa  <- 1
                   return(c(c(cc1[ind1[aa2]],cc1[ind2[aa2]]),0))
   				  } 
            })
            cor2EX <- TRUE
           }
          }
          ## pvalues second matrix
          cor3EX <- FALSE
          if(class(TRI2r)=="matrix"){
             if(dim(TRI2r)[1]>0){
                cor3 <- apply(as.matrix(1:dim(TRI2r)[1]),1,function(PL){
                  print(PL)
                  cc1 <- TRI2r[PL,]
                  aa <- try(abs(cov2cor(solve(cov2cor(cor(D2[,cc1]))))))
                  if(class(aa)!="try-error"){
                   aa2 <- which.min(lowerTri(aa))
                   aa  <- min(aa)
                   return(c(c(cc1[ind1[aa2]],cc1[ind2[aa2]]),pnorm(ztransf(aa)*sqrt(N-5))))
   				  }			  
                  else
   				  {
   				   aa2 <- 1
                   aa  <- 1
                   return(c(c(cc1[ind1[aa2]],cc1[ind2[aa2]]),0))
   				  }
            })
            cor3EX <- TRUE
           }
          }
          
          ## New adjacency matrices
          if(cor2EX){
           edgesToConsider[[1]]        <- cbind(cor2[1,],cor2[2,])
           pvaluesToConsider[[1]]      <- 1-cor2[3,]
           whzero                      <- c(which(pvaluesToConsider[[1]]>(alpha)))
           edgesToZero1                <- edgesToConsider[[1]][whzero,]
           if(!is.matrix(edgesToZero1)) edgesToZero1 <- t(as.matrix(edgesToZero1))
           A1.new[edgesToZero1]        <- 0
           A1.new[edgesToZero1[,2:1]]  <- 0
           if(equalBetweenCond){
            edgesToZero3 <- edgesToZero1[A1[edgesToZero1]==1& A2[edgesToZero1]==1,]
            if(class(edgesToZero3)!="matrix")
             edgesToZero3 <- t(as.matrix(edgesToZero3))
             A2.new[edgesToZero3]        <- 0
             A2.new[edgesToZero3[,2:1]]  <- 0
           }
          }
          if(cor3EX){
           edgesToConsider[[2]]        <- cbind(cor3[1,],cor3[2,])
           pvaluesToConsider[[2]]      <- 1-cor3[3,]
           whzero2                     <- which(pvaluesToConsider[[2]]>(alpha))
           edgesToZero2                <- edgesToConsider[[2]][whzero2,]
           if(!is.matrix(edgesToZero2)) edgesToZero2 <- t(as.matrix(edgesToZero2))
           A2.new[edgesToZero2]        <- 0
           A2.new[edgesToZero2[,2:1]]  <- 0
           if(equalBetweenCond){
            edgesToZero3 <- edgesToZero2[A1[edgesToZero2]==1& A2[edgesToZero2]==1,]
            if(class(edgesToZero3)!="matrix")
             edgesToZero3 <- t(as.matrix(edgesToZero3))
            A1.new[edgesToZero3]        <- 0
            A1.new[edgesToZero3[,2:1]]  <- 0
           }
          }
          
          return(list(A1.new=A1.new,A2.new=A2.new,edgesToConsider=edgesToConsider,
                      pvaluesToConsider=pvaluesToConsider,TRI1r=TRI1r,TRI2r=TRI2r))
}

asMatrix2 <- function(M){
 if(!is.matrix(M))
  {
   if(length(M)>0){
     M <- t(as.matrix(M))
   }
  }  
  return(M)
}
