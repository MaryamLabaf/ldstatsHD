amseLambdaSelection <- function(obj, pathIni, y, generator = c("subsampling", "montecarlo"),
								pB = 0.7, nite = 10, method = "mb", from = 1, until = NULL, 
								distF = c("correlation","shortPath"), oneByone = FALSE, 
                           		many = 3)
{
  ## checks 
  if (class(obj) != "huge" & class(obj) != "tiger")
   stop("obj must be a huge or tiger object")
  if(missing(y))
       stop("y argument (with dataset used) needs to be specified for amseLambdaSelection")
  
  if(missing(pathIni))
       stop("pathIni argument (with initial path) needs to be specified for amseLambdaSelection")
  else{
   if(dim(pathIni)[1]!=dim(obj$path[[1]])[1])
       stop("pathIni dimension is not equal to object path dimensions")
  }  	 

  if(!is.null(until)){
   if(until < 2 | until > length(obj$path))
	stop("UNTIL must be between 2 and the length of used lambdas for estimation")
  }
  else
   until <- length(obj$path)
  if(from < 0 | from >= until)
   stop("from must be larger than zero and smaller than until")
  
  pdistF <- c("correlation", "shortPath", "degrees", "corr01")
  if(length(distF) > 1)
   warning("distF attribute length is larger than one. Only the first component will be used")
  distF <- pdistF[pmatch(distF[1], pdistF)]
  if (is.na(distF)) stop("distF is not well define. It must be selected from \"correlation\", \"shortPath\", \"corr01\" or \"degrees\" ")
 
  pmethod <- c("mb", "glasso", "tiger")
  if(length(method) > 1)
   warning("method attribute length is larger than one. Only the first component will be used")
  method <- pmethod[pmatch(method[1], pmethod)]
  if (is.na(method)) stop("method is not well define. It must be selected from \"mb\",  \"glasso\" or \"tiger\" ")   
   
  pgenerator <- c("subsampling", "montecarlo")
  if(length(generator) > 1)
   warning("generator attribute length is larger than one. Only the first component will be used")
  generator <- pgenerator[pmatch(generator[1], pgenerator)]
  if (is.na(generator)) stop("generator is not well define. It must be selected from \"subsampling\" or \"montecarlo\" ")   

  if (pB < 0 | pB > 1)
   stop("pB is not well defined. Must be between 0 and 1")
	
  if (nite <= 0)
   stop("nite is not well defined. Must be larger than 0")
	
  if(oneByone & many > (until - from + 1))
   stop("many is not well defined. Must be smaller than until-from + 2")

 
  ## Initialization 
  PATHS   	<- obj$path
  LAMBDA  	<- obj$lambda
  P       	<- dim(PATHS[[1]])[1]
  VAR       <- rep(0,length(LAMBDA))
  BIASS     <- rep(0,length(LAMBDA))
  BS        <- rep(0,length(LAMBDA))
  N			<- 	dim(y)[1]  
  B         <- floor(N * pB)

  if(distF == "correlation") dist2     <- graphCorr(pathIni)
  if(distF == "shortPath")   dist2     <- graphDist(pathIni)
  if(distF == "degrees")     dist2     <- degrees(pathIni)
  if(distF == "corr01")      dist2     <- ifelse(graphCorr(pathIni) == 1, 0, 1)
                
   if(generator == "montecarlo")
   {
   	 D <- rep(0,P)
     COEF <- apply(as.matrix(1:P),1,function(j){
     if(sum(pathIni[j,] != 0)>0){
     	lm1 			   <- lm(y[,j] ~ y[,pathIni[j,]!=0]-1)
		D[pathIni[j,]!=0] <- lm1$coeff
		DIAG <- 1/var(lm1$resid)
	 }
   	 else DIAG <- 1
	 D   	<- -D * DIAG
     D[j]	<- DIAG
	return(D)
   })
   T1 		   <- forceSymmetric(COEF,"corMatrix")
   VV 		   <- eigen(as.matrix(T1))
   M  		   <- N
   DK          <- M/(2) * (-VV$val + sqrt(VV$val^2 + 4 * 1/M))
   AA          <- t(apply(VV$vec,1,function(x) x*DK))
   T2          <- AA %*% t(VV$vec)
   MAT 		   <- cov2cor(T2)
  }
  
  pb        <- txtProgressBar(min = 0, max = nite, style = 3)

  cont    <- 1
  PATHSd  <- vector(length(LAMBDA[from:length(LAMBDA)]),mode="list")
      
  ## AMSE
  if(generator == "subsampling"){
	  while(cont <= nite){
		  setTxtProgressBar(pb, cont)
		  IND    <- sample(1:dim(y)[1],B)
		  if (method == "mb" || method == "glasso")
		  {
			  if(oneByone){
				SEQ <- seq(from, until, by = many)
				for( i in SEQ[-length(SEQ)]){
				 print(i)
				 PATHSs <- huge(y[IND,], method = method, lambda=LAMBDA[i:min(until,(i+many-1))])$path
				 for(j in 1:many){
				   if(sum(PATHSs[[j]])>2)
				   {
					if(distF == "correlation") dist1     <- graphCorr(PATHSs[[j]])
					if(distF == "shortPath")   dist1     <- graphDist(PATHSs[[j]])
					if(distF == "degrees")     dist1     <- degrees(PATHSs[[j]])
					if(distF == "corr01")      dist1     <- ifelse(graphCorr(PATHSs[[j]]) == 1,0,1)
				
					if(cont == 1)  PATHSd[[i+j-1]]  <- sum((dist1-dist2)^2)
					if(cont > 1)   PATHSd[[i+j-1]]  <- PATHSd[[i+j-1]] + sum((dist1-dist2)^2)
				   }
				   else
					PATHSd[[(i+j-1)]]  <- 9999
				 }
			   }
			 }else
			   PATHSs <- huge(y[IND,], method = method, lambda = LAMBDA[from:until])$path
		   }
	   
		   if (method == "tiger")
			  PATHSs <- camel.tiger(y[IND,], method = "slasso", lambda = LAMBDA[from:length(LAMBDA)])$path
	   
		  if(!oneByone){
		   for(i in 1:(until-from+1))
		   {
				 if(sum(PATHSs[[i]]) > 2)
				   {
					if(distF == "correlation") dist1     <- graphCorr(PATHSs[[i]])
					if(distF == "shortPath")   dist1     <- graphDist(PATHSs[[i]])
					if(distF == "degrees")     dist1     <- degrees(PATHSs[[i]])
					if(distF == "corr01")      dist1     <- ifelse(graphCorr(PATHSs[[i]]) == 1,0,1)
					if(cont == 1)  PATHSd[[i]]  <- sum((dist1-dist2)^2)
					if(cont > 1)   PATHSd[[i]]  <- PATHSd[[i]] + sum((dist1-dist2)^2)
				   }
				   else
					PATHSd[[i]]    <- 99999
		   }
		   }
		   cont <- cont+1
	   }
	   AAG.COEF <- do.call(c,lapply(PATHSd, sum))/nite
   }
   if(generator == "montecarlo"){
      AAG.COEF <- 0
      while(cont <= nite){
   		  setTxtProgressBar(pb, cont)
		  y2 			<- mvrnorm( N, rep(0,P), MAT)
		  corrr 		<- cor(y2)
		  out4          <- huge(scale(y2), method = method, lambda=LAMBDA[from:until])
		  PATHSs        <- out4$path
		  DIST2 		<- list(length(LAMBDA[from:until]))
		  for(i in 1:length(LAMBDA[from:until])){
			if(distF == "correlation") DIST2[[i]]     <- graphCorr(PATHSs[[i]])
			if(distF == "shortPath")   DIST2[[i]]     <- graphDist(PATHSs[[i]])
			if(distF == "degrees")     DIST2[[i]]     <- degrees(PATHSs[[i]])
			if(distF == "corr01")      DIST2[[i]]     <- ifelse(graphCorr(PATHSs[[i]]) == 1,0,1)
		  }
		  AAG.COEF 		<- AAG.COEF + (unlist(lapply(DIST2, function(X) mean((X-dist2)^2))))
 		  cont 			<- cont + 1
     }
   }
   close(pb)
   
   ret.list   					<- list(opt.lambda =  LAMBDA[which.min(AAG.COEF) + from - 1], crit.coef = AAG.COEF, criterion = "A-MSE")
   attr(ret.list, "generator") 	<- generator
   ret.list$lambda    			<- obj$lambda[from:until]
   attr(ret.list, "bestpath")   <- obj$path[[ which(obj$lambda == ret.list[[1]][1]) ]]

   class(ret.list)  			<- "lambdaSelection"
   return(ret.list)    
}




 