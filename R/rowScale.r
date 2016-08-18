rowVars <- function(mat)  apply(mat,1,var)

rowScale<-function(mat){
  mat.x<-mat[rowVars(mat)>0,]
  new.mat<-t(scale(t(mat.x)))
  return(new.mat)
}
