vulLambdaSelection <- function(obj, loo = FALSE, subOut =  10, nite = 50)
{
        ## checks
	    if(class(obj) != "huge" & class(obj) != "tiger" & class(obj) != "wfgl")
 		 stop("object has to be of class huge, tiger or wfgl")
	    if(class(obj) == "wfgl")
		 obj$path <- lapply(obj$path,function(x) (x[[1]]|x[[2]])*1)
	    if(class(obj) == "wfgl")
		 obj$lambda <- obj$lambda1
		if(length(obj$lambda) < 5)
         stop("lambda sequence has to be of length 5 or higher")
		if(subOut <= 0 | subOut >= dim(obj$path[[1]])[1])
         stop("subOut has to be larger than 0 and smaller than p")
        if(nite < 2)
         stop("nite has to be larger than 2")         
		
        ## initialization 
        PATHS  		<- obj$path
        LAMBDA 		<- obj$lambda
        P 			<- dim(PATHS[[1]])[1]
		minDegree 	<- 2
		
		## vulnerability
        pb <- txtProgressBar(min = 0, max = length(LAMBDA), style = 3)
        HarmonicMeans <- apply(as.matrix(1:length(LAMBDA)), 1, function(k)
        {
           PATH     <- PATHS[[k]]
           setTxtProgressBar(pb, k)
           NODES.DEGREE    <- degrees(PATH)
           if(sum(NODES.DEGREE >= minDegree)> subOut)
           {
            E <- harmMeanDist(PATH, NODES.DEGREE)
            
            if (!loo){
                Es <- apply(as.matrix(1:nite),1,function(j)
                {
                   IND <- sample(1:P,floor(subOut))
                    harmMeanDist((PATH[,-IND])[-IND,], NODES.DEGREE[-IND])
                }) 
                Su <- sum((E-mean(Es))/E)
            }
            if (loo){
               Es <- apply(as.matrix(which(NODES.DEGREE >= 0)),1,function(j)
               harmMeanDist((PATH[,-j])[-j,],NODES.DEGREE[-j]))
               Su <- sum((E-Es)/E)
            }
            return(Su)
           }
           else 
              return(0)
         })
         
        lambdaVUL <- LAMBDA[which.min(HarmonicMeans)]
        close(pb)
        
        ret.list   		<- list(opt.lambda = lambdaVUL, crit.coef = HarmonicMeans, criterion = "VUL")
        ret.list$lambda <- LAMBDA
        
        attr(ret.list, "bestpath")   <- obj$path[[ which(obj$lambda == ret.list[[1]][1]) ]]
        
        class(ret.list) <- "lambdaSelection"

		return(ret.list)
}

