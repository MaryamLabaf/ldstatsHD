plot.eqCorTestByRows <- function(x, mains = c("AS CI", "max CI"), xlabs = c("",""), ylabs = c("",""), pch = "-", ownCols =TRUE, ...){
  if (!inherits(x, "eqCorTestByRows")) stop("use only with \"eqCorTestByRows\" object")
  testStatistic <- x$testStatistic
  if(ownCols){
   if(any(testStatistic == "AS"))
   {
    col1 <- rep(1, length(x$ciAS[,2]))
    col1[x$ciAS[,1]>0] <- 3
   }
   
   if(any(testStatistic == "max"))
   {
    col2 <- rep(1, length(x$ciMax[,2]))
    col2[x$ciMax[,1]>0] <- 3
   }
   
   if(length(testStatistic) == 2) par(mfrow=c(2,1))
   if(any(testStatistic == "AS"))
   {
    plot(x$ciAS[,2], type = "p", ylim=c(min(x$ciAS),max(x$ciAS)), pch = pch, xlab = xlabs[1], main = mains[1], ylab= ylabs[1], col = col1, ...)
    points(x$ciAS[,1], pch="--", col = col1, ...)
    for(k in 1:dim(x$ciAS)[1]) lines(c(k,k),c(x$ciAS[k,1],x$ciAS[k,2]), col = col1[k], ...)
    abline(h=0,col=2,lty=2)
   }
   if(any(testStatistic == "max"))
   {
    plot(x$ciMax[,2],type = "p", ylim=c(min(x$ciMax),max(x$ciMax)), pch = pch, xlab = xlabs[2], main = mains[2], ylab= ylabs[2], col = col2, ...)
    points(x$ciMax[,1], pch="--", col = col2, ...)
    for(k in 1:dim(x$ciMax)[1]) lines(c(k,k),c(x$ciMax[k,1],x$ciMax[k,2]), col = col2[k], ...)     
    abline(h = 0,col = 2, lty = 2)
   }
  }
  else{
   if(any(testStatistic == "AS"))
   {
    if(length(testStatistic) == 2) par(mfrow=c(2,1))
    plot(x$ciAS[,2], type = "p", ylim=c(min(x$ciAS),max(x$ciAS)), pch = pch, xlab = xlabs[1], main = mains[1], ylab= ylabs[1], ...)
    points(x$ciAS[,1], pch="--",...)
    for(k in 1:dim(x$ciAS)[1]) lines(c(k,k),c(x$ciAS[k,1],x$ciAS[k,2]),...)
    abline(h=0,col=2,lty=2)
   }
   if(any(testStatistic == "max"))
   {
    plot(x$ciMax[,2],type = "p", ylim=c(min(x$ciMax),max(x$ciMax)), pch = pch, xlab = xlabs[2], main = mains[2], ylab= ylabs[2], ...)
    points(x$ciMax[,1], pch="--",...)
    for(k in 1:dim(x$ciMax)[1]) lines(c(k,k),c(x$ciMax[k,1],x$ciMax[k,2]),...)     
    abline(h = 0, col = 2, lty = 2)
   }
  }
}