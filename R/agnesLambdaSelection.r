### Function agnesLambdaSelection selects the regularisation parameter lambda by maximising the AGNES coefficient.
### way defines if we want to use a variable subset approach (rand.sampling, int.sampling) or not (direct)
agnesLambdaSelection <- function(obj, way = "direct", nite = 10, subsvec = NULL,
                                eps = 0.05, until = NULL, minNodes = 30, 
                                distF = c("correlation","shortPath"))
{
      ## checks
	  if(class(obj) != "huge" & class(obj) != "tiger" & class(obj) != "wfgl")
 		stop("object has to be of class huge, tiger or wfgl")
	  if(class(obj) == "wfgl")
		obj$path <- lapply(obj$path,function(x) (x[[1]]|x[[2]])*1)
	  if(class(obj) == "wfgl")
		obj$lambda <- obj$lambda1
		
      ALFA <- 1 ## internal
      if(!is.null(until)){
       if(until < 2 | until > length(obj$path))
        stop("until must be between 2 and the length of used lambdas for estimation")
      }
      else
       until <- length(obj$path)
      
      pway <- c("direct", "rand.sampling", "int.sampling")
      if(length(way) > 1)
       warning("way attribute length is larger than one. Only the first component will be used")
   	  way <- pway[pmatch(way[1], pway)]
      if (is.na(way)) stop("way is not well define. It must be selected from \"direct\", \"rand.sampling\" or \"int.sampling\" ")
       
      pdistF <- c("correlation", "shortPath")
      if(length(distF) > 1)
       warning("distF attribute length is larger than one. Only the first component will be used")
   	  distF <- pdistF[pmatch(distF[1], pdistF)]
      if (is.na(distF)) stop("distF is not well define. It must be selected from \"correlation\" or \"shortPath\" ")
 
 	  if (length(nite) == 1 & way != "direct") nite <- rep(nite, until)
 	  if(way != "direct" & any(nite < 1))
  	    stop("nite elements must be at least 1") 
 	  
 	  if (is.null(subsvec) & way != "direct") subsvec <- rep(minNodes, until)
 	  if(way != "direct" & any(subsvec < 5))
  	    stop("subsvec elements must be at least 5") 
 	  
 	  if(way != "direct" & length(subsvec) < until)
  	    stop("length(subsvec) must be at least until") 
 	 
 	  if(eps <= 0)
 	   stop("eps must be larger than zero")
 
      ## Initialization
      PATHS   		<- obj$path
      LAMBDA  		<- obj$lambda
      P       		<- dim(PATHS[[1]])[1]
      MINT    		<- 300   #internal control
      
      AGL.COEF 		<- rep(0,length(LAMBDA))
      P        		<- dim(PATHS[[1]])[1]
      
      pb       <- txtProgressBar(min = 0, max = (until), style = 3)
      
      ## AGNES
      for (SI in 1:until)
      {
	    setTxtProgressBar(pb, SI)
        A         <- PATHS[[SI]]
        diag(A)   <- rep(1,dim(A)[1])

        ## Pearson correlation
        KS 	  <- degrees(A) 	
        IND   <- which(KS > 1)

        if (length(IND) > minNodes)
        {
           A        <- A[IND,IND]
           P2       <- length(IND)
           KS       <- KS[IND]
           CV.T     <- sd(KS)/mean(KS)
           if(distF == "correlation") DIST.MA  <- graphCorr(A,KS)
           if(distF == "shortPath")   DIST.MA  <- graphDist(A)
           
           if(way != "direct" & any(subsvec > P2))
           {
            warning("subsvec have elements larger than the number of effective nodes (which will be used instead")
            subsvec[subsvec > P2] <- P2
           }
           
           if (way == "rand.sampling")
           {
                DIST.AUX      	<- array(0,dim=c(P2,P2))
                suppressMessages(DIST.AUX[lower.tri(DIST.AUX)] <- DIST.MA)
                DIST.AUX      	<- DIST.AUX +t(DIST.AUX)
                ALFA.TEMP 		<- 0; AGL.COEF.TEMP <- 0
                nite2     		<- nite[SI]
                J         		<- 1
                while (J <= nite2)
                {
                        IND2          <- sample(seq_len(P2), min(P2,subsvec[SI]))
                        SORT          <- agnes(DIST.AUX[IND2,IND2], diss = TRUE)
                        AUX           <- DIST.AUX[IND2,IND2]
                        AGL.COEF.TEMP <- AGL.COEF.TEMP +  SORT$ac
                        ALFA.TEMP     <- ALFA.TEMP + length(IND)/P2
                        J             <- J+1
                }
                AGL.COEF[SI] <- AGL.COEF.TEMP / nite2 * (1- ALFA * ALFA.TEMP/nite2 * (1- P2/P))
           }
           if (way == "int.sampling")
           {
                AGL.COEF.TEMP 	<- 0
                ALFA.TEMP     	<- 0
                nite2         	<- nite[SI] 
                J  				<- 1
                KA				<- 1
                DIST.AUX      	<- array(0,dim=c(P2,P2))
                suppressMessages(DIST.AUX[lower.tri(DIST.AUX)] <- DIST.MA)
                DIST.AUX      	<- DIST.AUX +t(DIST.AUX)
                ALFA.TEMP 		<- 0; AGL.COEF.TEMP <- 0
                J         		<- 1
                     
                while (J <= nite2 & KA < 100 )
                {
					  IND2          <- sample(seq_len(P2),min(P2,subsvec[SI]))
					  IND2          <- apply(as.matrix(IND2),1,function(j) which(A[j,]!=0))
					  if (class(IND2)=="list") IND2 <- unique(do.call(c,IND2))
						IND2 <- sort(sample(IND2,min(MINT,length(IND2))))
						CV.E <- sd(KS[IND2])/mean(KS[IND2])
						if (abs(CV.E-CV.T)/CV.T < eps){
						  SORT          <- agnes(DIST.AUX[IND2,IND2],diss=TRUE)
						  AUX           <- DIST.AUX[IND2,IND2]
						  AGL.COEF.TEMP <- AGL.COEF.TEMP +  SORT$ac
						  ALFA.TEMP     <- ALFA.TEMP + length(IND)/P2
						  J             <- J+1
						}          
						KA <- KA+1                           
				}                
                AGL.COEF[SI] <- AGL.COEF.TEMP / nite2 * (1- ALFA * ALFA.TEMP/nite2 * (1- P2/P))
            }
            if (way == "direct"){
                SORT          <- agnes(DIST.MA, diss=TRUE)
                AGL.COEF[SI]  <- SORT$ac * (1-ALFA*(1- P2/P))
            }
         }
         else AGL.COEF[SI] <-0
         setTxtProgressBar(pb, SI)
        }
        close(pb)

        ret.list   		 <- list(opt.lambda = LAMBDA[which.max(AGL.COEF)], crit.coef = AGL.COEF, criterion = "AGNES")
        ret.list$lambda  <- obj$lambda
        attr(ret.list, "bestpath") <- obj$path[[ which(obj$lambda == ret.list[[1]][1]) ]]
  
        class(ret.list)  <- "lambdaSelection"

		return(ret.list)
}


