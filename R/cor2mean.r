### function that calculates the average unadjusted squared correlations of the each row
###of the correlation matrix including the diagonal
cor2mean <- function(mat){
  n <- dim(mat)[2]
  p <- dim(mat)[1]
  cor2.vec <- rep(0,p)
  mat.r <- rowScale(mat)
  p <- dim(mat.r)[1]
  A <- mat.r%*%(t(mat.r) %*%(mat.r))

  for (i in 1:p) cor2.vec[i] <- 1/(p*(n-1)^2)* A[i,]%*%mat.r[i,]

  return(cor2.vec)
}
