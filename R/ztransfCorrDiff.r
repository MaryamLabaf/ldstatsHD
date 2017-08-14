ztransfCorrDiff  <- function(D1, D2, thetaKnown=NULL, biasCorr=NULL ){
    N     <- dim(D1)[1]
    P     <- dim(D1)[2]
    ID    <- 1:N
#    CJ    <- cor(cbind(D1,D2))
#    C1   <- lowerTri(CJ[1:P,1:P])
#    C2   <- lowerTri(CJ[1:P+P,1:P+P])
	C1   <- lowerTri(cor(D1))
	C2   <- lowerTri(cor(D2))
    if(is.null(thetaKnown)){
      Ca	<- cor(D1,D2)	 	
      CC    <- diag(Ca)
      C12d1 <- do.call(c,apply(as.matrix(1:P),1,function(i) rep(CC[i],P-i)))
      C12d2 <- do.call(c,apply(as.matrix(1:(P-1)),1,function(i) CC[(i+1):P]))
#      C12   <- lowerTri(CJ[1:P,1:P+P])
#      C21   <- lowerTri(CJ[1:P+P,1:P])
	   C12   <-lowerTri(Ca)
	   C21	<- lowerTri(t(Ca))	 	
      THETA <- pmin(0.95,covCrosCor((C1+C2)/2, (C1+C2)/2, C12, C12d1, C12d2, C21))
    }
    else
     THETA <- thetaKnown

    if(!(is.null(biasCorr))) THETA <- THETA + biasCorr    
    ZtransfD <- ztransf(C1)*sqrt(N-3) - ztransf(C2)*sqrt(N-3)
    ListRET  <- list(Ts = ZtransfD/sqrt(2-2*THETA), theta=THETA, Ts2= ZtransfD)
    return(ListRET)
}

thetaFromCor  <- function(COR, P){
      C1   <- lowerTri(COR[1:P,1:P])
      C2   <- lowerTri(COR[1:P+P,1:P+P])
      CC    <- diag(COR[1:P,1:P+P])
      C12d1 <- do.call(c,apply(as.matrix(1:P),1,function(i) rep(CC[i],P-i)))
      C12d2 <- do.call(c,apply(as.matrix(1:(P-1)),1,function(i) CC[(i+1):P]))
      C12   <- lowerTri(COR[1:P,1:P+P])
      C21   <- lowerTri(COR[1:P+P,1:P])
      THETA <- pmin(0.95,covCrosCor(C1, C2, C21, C12d1, C12d2, C12))
   
    return(THETA)
}

ztransfCorrDiffsub  <- function(D1, D2, row, thetaKnown = NULL,  diagCor ){
    N     <- dim(D1)[1]
    P     <- dim(D1)[2]
    ID    <- 1:N
    ulr   <- 1:length(row)
    C1    <- as.numeric(t(cor(D1[,row],D1[,])))[- ((ulr-1)*P +row)]
    C2    <- as.numeric(t(cor(D2[,row],D2[,])))[- ((ulr-1)*P +row)]
    if(is.null(thetaKnown)){
      C12   <- as.numeric(t(cor(D1[,row],D2[,])))[- ((ulr-1)*P +row)]
      C21   <- as.numeric(t(cor(D2[,row],D1[,])))[- ((ulr-1)*P +row)]
      C12d1 <- apply(as.matrix(row),1,function(j) rep(cor(D1[,j],D2[,j]),P))[- ((ulr-1)*P +row)]
      C12d2 <- rep(diagCor,length(row))[- ((ulr-1)*P +row)]#[-row]
      THETA <- pmin(0.95,covCrosCor((C1+C2)/2, (C1+C2)/2, C12, C12d1, C12d2, C21))
    }
    else
     THETA <- thetaKnown

    ZtransfD <- ztransf(C1)*sqrt(N-3) - ztransf(C2)*sqrt(N-3)
    ListRET  <- list(Ts = ZtransfD/sqrt(2-2*THETA), theta=THETA, Ts2= ZtransfD)
    return(ListRET)
}