pcLambdaSelection <- function(obj)
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

        ## initialization
        PATHS  		<- obj$path
        LAMBDA 		<- obj$lambda
        P 			<- dim(PATHS[[1]])[1]
        minNodes	<- 3
        
        ## path connectivity
        pb      	<- txtProgressBar(min = 0, max = length(LAMBDA), style = 3)
        SUMdist 	<- try(apply(as.matrix(1:length(LAMBDA)),1, function(k)
        {
                PATH            <- PATHS[[k]]
                setTxtProgressBar(pb, k)
                NODES.DEGREE    <- degrees(PATH)
                if (any(NODES.DEGREE >= minNodes ))
                {
                  NODESINGRAPH  <- which(NODES.DEGREE > 0)
                  g             <- graph.adjacency(Matrix(PATH[NODESINGRAPH,NODESINGRAPH], sparse = TRUE), mode = "undirected")
                  g             <- igraph.to.graphNEL(g)
                  DIST          <- as.dist(johnson.all.pairs.sp(g))
                  WHICH.DIST    <- sum(DIST[is.finite(DIST)])
                  return(WHICH.DIST)
                }else
                  return(0)
        } ))
                 
        LENG     <- sum(SUMdist > 0)
        SUMdist  <- SUMdist[1:LENG]
        DSUMdist <- diff(SUMdist)

        runAverage   <- apply(as.matrix(1:(LENG-1)),1,function(i){
                          mean(DSUMdist[i:(LENG-1)])
                         })
        lambdaPC 	 <- LAMBDA[which.max(abs(DSUMdist[1:(LENG-1)]/runAverage))+1]
        close(pb)
      
        ret.list   		<- list(opt.lambda = lambdaPC, crit.coef = SUMdist, criterion = "PC")
        ret.list$lambda <- LAMBDA[1:LENG]
        
        attr(ret.list, "bestpath") <- obj$path[[ which(obj$lambda == ret.list[[1]][1]) ]]

        class(ret.list) <- "lambdaSelection"

		return(ret.list)
}
