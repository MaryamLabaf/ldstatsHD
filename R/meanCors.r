
meanCors <- function(mat, scaleT = TRUE, sampleFrom = FALSE, nitSamp = 1000){
  n <- dim(mat)[2]
  p <- dim(mat)[1]
  rep.v <- rep(1,p)
  if(scaleT)  mat.r <-rowScale(mat)
  else  mat.r <- mat
  A2 <-(t(mat.r) %*% t(t(rep.v)))
  A3 <- t(rep.v) %*% mat.r
  cor1.vec <- 1/(p^2*(n-1)) * A3%*%A2

  if(scaleT){
    mat.r <-rowScale(mat^2)
    A2 <-(t(mat.r) %*% t(t(rep.v)))
    A3 <- t(rep.v) %*% mat.r
    cor2.vec <- (1/(p^2*(n-1)) * A3%*%A2)
  }
  else{
    if(sampleFrom) mat   <- mat[sample(1:p,nitSamp),]
    mat.r <- mat^2
    mat.r <- (mat.r-apply( mat.r,1,mean))

    p<-dim(mat.r)[1]
    rep.v <-rep(1,p)
    A2 <-(t(mat.r) %*%t(t(rep.v)))
    A3 <- t(rep.v) %*% mat.r
    cor2.vec <- (1/(p^2*(n-1)) * A3%*%A2)/2
  }
  rm(A2)
  return(c(cor1.vec, cor2.vec))
}
