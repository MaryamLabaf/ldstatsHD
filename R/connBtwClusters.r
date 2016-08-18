
######### Random perturbation between clusters in the partial correlation ######
connBtwClusters  <- function(PREC.MAT, perturb.clust=0)
{

  PREC.MAT            <- as.matrix(PREC.MAT)
  PREC.MAT.AUX        <- PREC.MAT
  diag(PREC.MAT.AUX)  <- 0
  sumNotZero          <- sum(PREC.MAT.AUX!=0)/2

  diag(PREC.MAT.AUX)  <- diag(PREC.MAT)
  if(perturb.clust > 0)
  {
    WHI1 <- which(PREC.MAT == 0, arr.ind=TRUE)
    WHI1 <- WHI1[WHI1[,1]<WHI1[,2],]
    IND1 <- sample(1:dim(WHI1)[1], round(sumNotZero * perturb.clust))

    SignInEdges     <- rbinom(length(IND1),1,0.5)
    SignInEdges     <- ifelse(SignInEdges==0,-1,1)
    IntenInEdges    <- runif(length(IND1), min(abs(PREC.MAT[PREC.MAT!=0])),
				sort(abs(suppressMessages(PREC.MAT[lower.tri(PREC.MAT)])),decreasing = TRUE)[5]) * SignInEdges
    
    PREC.MAT.AUX[WHI1[IND1,]]    <- IntenInEdges
    PREC.MAT.AUX[WHI1[IND1,2:1]] <- IntenInEdges
    
  }
   delta  <- 0
   eigval <- eigen(as.matrix(PREC.MAT.AUX))$val

   while(max(eigval)/min(eigval) > dim(PREC.MAT)[1] | min(eigval) <0)
   {
    delta        <- delta + 0.01
    PREC.MAT.AUX <- PREC.MAT.AUX + delta * diag(rep(1,dim(PREC.MAT.AUX)[1]))
    eigval <- eigen(as.matrix(PREC.MAT.AUX))$val
  }
  return(PREC.MAT.AUX)

}


