degrees <- function(A){
 P   <- dim(A)[1]
 deg <- numeric(P)
 for(i in 1:P) deg[i] <- sum(A[i,])
 return(deg)
}