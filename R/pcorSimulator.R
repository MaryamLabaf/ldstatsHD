
########################## Partial correlation and data simulator ##############
pcorSimulator  <- function(nobs, nclusters, nnodesxcluster, pattern = "powerLaw", 
                           low.strength = 0.5, sup.strength = 0.9, nhubs = 5, 
                           degree.hubs = 20, nOtherEdges = 30, alpha = 2.3, plus = 0, 
                           prob = 0.05, perturb.clust = 0, mu = 0,
                           probSign = 0.5, seed = sample(10000, nclusters))
{

  ###### Checks
  cPS <- controlsPcorSimulator(nobs, nclusters, nnodesxcluster, pattern, 
                      low.strength, sup.strength, nhubs, degree.hubs, nOtherEdges, 
                      alpha, plus, prob, perturb.clust, mu, probSign, seed)

  ## GRAPH CONSTRUCTION
  orderNodes1 <- c(0, cumsum(cPS$nnodesxcluster))
  MOD  <- lapply(as.matrix(1:cPS$nclusters),function(i)
          graphStructure(nhubs = cPS$nhubs[i], alpha = cPS$alpha, nnodes = cPS$nnodesxcluster[i], 
                         degree.hubs = cPS$degree.hubs[i], low.strength = cPS$low.strength, 
                         sup.strength = cPS$sup.strength, nOtherEdges = cPS$nOtherEdges[i],
                         orderNodes = orderNodes1[i], probSign = cPS$probSign, pattern = cPS$pattern, 
                         seed = cPS$seed[i], plus = cPS$plus, prob = cPS$prob[i]))

  EdgesToChange1 	<- do.call(rbind, lapply(MOD,function(x) x$edgesToChange))         
  hubs          	<- as.numeric(names(table(EdgesToChange1)))[which(table(EdgesToChange1) 
  										>= max(degree.hubs))]
  omega       		<- as.matrix(do.call(bdiag, lapply(MOD,function(x) x$omega)))
  
  ## Between clusters connections 
  if(cPS$perturb.clust > 0 & cPS$nclusters > 1){
    omega  <- connBtwClusters(omega, perturb.clust = cPS$perturb.clust)
  }
  
  
  ## Finding covariance matrix and path
  rm(list = ".Random.seed", envir = globalenv()) 
  covMat 		<- pseudoinverse(omega)
  
  path        	<- ((omega)!=0)*1
  diag(path)  	<- 0

  ## Generating observations and return
  y   <- mvrnorm(nobs, cPS$mu, covMat)  

  obj <- list(y = y, hubs = hubs, edgesInGraph = EdgesToChange1, omega = omega, 
              covMat = covMat, path = path, pattern = cPS$pattern)
  class(obj) <- "pcorSim"            

  return(obj)
}

