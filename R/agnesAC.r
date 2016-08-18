agnesCoef <- function(A)
{
#	if(ALFA < 0 | ALFA > 1) stop("ALFA must be between 0 and 1")
 	ALFA		  <- 1		
    P             <- dim(A)[1]
    KS            <- degrees(A)
    NODESINGRAPH  <- which(KS > 0)
    DIST.MA       <- rep(0, P * (P-1) / 2)
    if(length(NODESINGRAPH)>0)
    {
        IND     <- which(KS > -1)
        P2      <- length(IND)
        A       <- A[IND,IND]
        DIST.MA <- graphCorr(A)

        DIST.AUX      <- array(0, dim = c(P2, P2))
        suppressMessages(DIST.AUX[lower.tri(DIST.AUX)] <- DIST.MA)
        DIST.AUX      <- DIST.AUX + t(DIST.AUX)
        SORT          <- agnes(DIST.AUX)
    }
    return(SORT$ac * (1 - ALFA * (1 - P2/P)))
}
