#gene subindexes in a lower triangular matrix
lowerTriMatInd<- function(x, P){
    i1 <- ceiling(((P-1/2)-sqrt((P-1/2)^2 -2*x)))
    i2 <- P-(i1*P-i1*(i1+1)/2-x)
    return(list(i1,i2))
}

