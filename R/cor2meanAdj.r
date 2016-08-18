###function that calculates the average ADJUSTED squared correlations of each row
###of the correlation matrix EXCLUDING the diagonal

cor2mean.adj <- function(mat){
  
  n <- dim(mat)[2]
  p <- dim(mat)[1]
  
  without.diag <- p/(p-1)*cor2mean(mat)-1/(p-1)
  return((n-1)/(n-2)*without.diag-1/(n-2))
}
