saddlePointFunction <- function(Exc, NeR, N, u, ttmin=-1.45, ttmax=.4999,lengthTT=5000, kmin=0, kmax=100, excAdj = FALSE){
 tt      <- seq(ttmin, ttmax,length.out=lengthTT)
 dn 		<- dnorm(u*(1-2*tt)^(1/2))
 pn1 	<- pnorm(-u)
 pn2 	<- pnorm(-u*(1-2*tt)^(1/2))
 K0a  	<- -1/2*log(1-2*tt) + log(pn2)-log(pn1)
 K1a 	<- 1/(1-2*tt) + 1/(1-2*tt)^(1/2) * u * dn/pn2
 K2a 	<- 2/(1-2*tt)^2 + u/sqrt(1-2*tt) * dn/pn2 * (u^2 + 1/(1-2*tt) - u * dn/pn2 *1/sqrt(1-2*tt))

 if(length(Exc)==1) pval <- numeric(kmax+1)
 if(length(Exc)>1) pval <- array(0,dim=c(length(Exc),kmax+1))

 if(kmin==0){
      if(length(Exc)==1) pval[i,1] <-  dbinom(0, N, 2*pn1) * 1
      if(length(Exc)>1)  pval[,1]  <-  dbinom(0, N, 2*pn1) * 1
 }

 if(kmin<=1){ 
   if(length(Exc)==1)
    if(Exc>0)  pval[2] <-  dbinom(1, N, 2*pn1) * (1- (pnorm(-sqrt(Exc)))/pn1)
   if(length(Exc)>1)
    {
      for(i in which( Exc > u^2)){
       pval[i,2] <-  dbinom(1, N, 2*pn1) * (1- (pnorm(-sqrt(Exc[i])))/pn1)
    }  
  }    
 } 
 for(k in max(2,kmin):kmax){
   print(k)
   Nu <- k
   fNu <- ((Nu*2*pi*K2a))^(-1/2) * exp(Nu*K0a - tt*K1a*Nu)
   fNu <- fNu/max(fNu)
   FNu <- cumsum(diff(K1a)*fNu[-1])/sum(diff(K1a)*fNu[-1])

   if(length(Exc)==1)
   {
     Excmean <- Exc/Nu
     if(Excmean > u^2)
     {
       wh   <- which(Excmean < K1a)
       pval[k+1] <- dbinom(Nu, N, 2*pn1) * ifelse(length(wh)>=1, FNu[wh[1]],1)
    }
   }
   else{
    for(i in which( Exc > ((u^2)*2))){
     Excmean <- Exc[i]/Nu
     if(Excmean > u^2)
     {
       wh        <- which(Excmean < K1a)
       pval[i,k+1] <- dbinom(Nu, N, 2*pn1) * ifelse(length(wh)>=1, FNu[wh[1]],1)
    }
   }
  }
 }	
 return(pval)
}
