### Function graphCorr finds the dissimilarity matrix d= 1- \frac{\eta_{ij}}{\sqrt{k_ik_j}}
### where eta_ij are shared neighbours btw nodes i and j and k_i is the degree of the node i
graphCorr <- function(A, nodesDegree = NULL)
{
        P    <- dim(A)[1]
        N    <- P * (P-1)/2
        N1   <- rep(0,N)
        N2   <- rep(0,N)
        
	    if(is.null(nodesDegree))  nodesDegree  <- degrees(A)
	    
        SIMIL2den   			<- nodesDegree %*% t(nodesDegree)
        A           			<- as.matrix(A)
        SIMIL2      			<-  suppressMessages((tcrossprod(A))[lower.tri(A)])/sqrt((SIMIL2den)[lower.tri(A)])
        SIMIL2[is.nan(SIMIL2)] 	<- 0
        DISSIM 					<- 1 - SIMIL2
         
        rm(SIMIL2); rm(SIMIL2den)
        return(DISSIM)
}

