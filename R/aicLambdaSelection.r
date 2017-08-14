 
aicAndbicLambdaSelection <- function(obj, y,  criterion = c("AIC", "BIC", "eBIC"),phi=1)
{
        ## checks
        if(class(obj) != "huge" & class(obj) != "tiger")
         stop("object has to be of class huge or tiger")
        if(length(obj$lambda) < 5)
         stop("lambda sequence has to be of length 5 or higher")
	    if(missing(y))
         stop("y argument (with dataset used) needs to be specified for aicAndbicLambdaSelection")
	    if(is.null(obj$icov) | is.null(obj$cov))
         stop("obj$icov and obj$cov have to be defined for aicAndbicLambdaSelection")
      
        pcriterion <- c("AIC", "BIC", "eBIC")
        criterion <- pcriterion[pmatch(criterion[1], pcriterion)]
		if (is.na(criterion)) stop("criterion is not well define. It must be selected from \"AIC\" or \"BIC\" or  \"eBIC\" ")

        ## initialization
        LAMBDA  <- obj$lambda
        N       <- dim(y)[1]
        pb      <- txtProgressBar(min = 0, max = length(LAMBDA), style = 3)
        M		<- dim(y)[2]*(dim(y)[2]-1)
        ## AIC and BIC
        AICs <-apply(as.matrix(1:length(LAMBDA)), 1, function(k){
          setTxtProgressBar(pb, k)
          aa 	<- (y)%*%as.matrix(obj$icov[[k]])
          AA2 	<-  2*N*(0.5 * log(det(as.matrix(obj$cov[[k]])))) - 2*sum(aa*y)
          BB 	<- sum(unlist(lapply(AA2,mean)))
          AL 	<- sum(obj$path[[k]]) 
          c(AL + BB, AL/2*log(N) + BB, AL/2*log(N) + BB + 2*phi*lchoose(M,sum(obj$path[[k]])))
        })
        
        lambdaAIC  <- LAMBDA[which.min(AICs[1,])]
        lambdaBIC  <- LAMBDA[which.min(AICs[2,])]
        lambdaeBIC <- LAMBDA[which.min(AICs[3,])]
        close(pb)

        if(length(criterion) == 2)
        {
         ret.list    	<- list(opt.lambda = c(lambdaAIC, lambdaBIC), 
                            crit.coef = list(AIC.COEF = AICs[1,], BIC.COEF = AICs[2,]), 
                            criterion = c("AIC", "BIC"))
         attr(ret.list, "bestpath2") <- obj$path[[ which(obj$lambda == ret.list[[1]][2]) ]]
		 attr(ret.list, "bestpath") <- obj$path[[ which(obj$lambda == ret.list[[1]][1]) ]]
		}
		else
        {
         if(criterion == "AIC")
         {
           ret.list    	<- list(opt.lambda = c(lambdaAIC), 
                            crit.coef = AICs[1,], 
                            criterion = "AIC")
		 }
		 if(criterion == "BIC")
         {
           ret.list    	<- list(opt.lambda = c(lambdaBIC), 
                            crit.coef = AICs[2,],
                            criterion = "BIC")
		 }
		 if(criterion == "eBIC")
         {
           ret.list    	<- list(opt.lambda = c(lambdaeBIC), 
                            crit.coef = AICs[3,],
                            criterion = "eBIC")
		 }
		 attr(ret.list, "bestpath") <- obj$path[[ which(obj$lambda == ret.list[[1]][1]) ]]
        }
        ret.list$lambda  <- LAMBDA		
        class(ret.list)  <- "lambdaSelection"
		return(ret.list)
}

