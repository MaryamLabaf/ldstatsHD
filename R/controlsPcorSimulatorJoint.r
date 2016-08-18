controlsPcorSimulatorJoint <- function(nobs, nclusters, nnodesxcluster, pattern, 
                      diffType, dataDepend, low.strength, sup.strength, pdiff, nhubs, 
                      degree.hubs, nOtherEdges, alpha, plus, prob, perturb.clust, mu, 
                      diagCCtype, diagNZ.strength, mixProb, probSign, seed, exactZeroTh)
{
  ## prospective patterns
  ppaterns <- c("hubs", "powerLaw", "random")
  if(length(pattern)>1) 
   warning("pattern attribute length is larger than one. Only the first component will be used")
  pattern <- ppaterns[pmatch(pattern[1], ppaterns)]
  if (is.na(pattern)) stop("pattern is not well define. It must be selected from \"hubs\", \"powerLaw\" or \"random\" ")
    
  ## prospective diffType
  pdiffType <- c("cluster", "random", "mixed")
  if(length(diffType)>1) 
   warning("diffType attribute length is larger than one. Only the first component will be used")
  diffType <- pdiffType[pmatch(diffType[1], pdiffType)]
  if (is.na(diffType)) stop("diffType is not well define. It must be selected from \"cluster\", \"random\" or \"mixed\" ")
    
  ## Data dependence
  pdataDepend <- c("ind", "add", "mult", "diagOmega")
  if(length(dataDepend)>1) 
   warning("dataDepend attribute length is larger than one. Only the first component will be used")
  dataDepend <- pdataDepend[pmatch(dataDepend[1], pdataDepend)]
  if (is.na(dataDepend)) stop("dataDepend is not well define. It must be selected from \"ind\", \"add\", \"mult\" or \"diagOmega\" ")
   
   ## nobs, nclusters and nnodesxclustern  
  if(nobs < 10) stop("nobs must be larger 10")
  if(nclusters < 1) stop("nclusters is not well defined. At least one cluster is needed") 
  if(length(nnodesxcluster)!= nclusters) 
    stop("nnodesxcluster must be a vector of length defined by nclusters")
  
  ## Data dependence: diagonal
  pdiagCCtype <- c("dicot", "beta13")
  if(length(diagCCtype)>1) 
   warning("diagCCtype attribute length is larger than one. Only the first component will be used")
  diagCCtype <- pdiagCCtype[pmatch(diagCCtype[1], pdiagCCtype)]
  if (is.na(diagCCtype)) stop("diagCCtype is not well define. It must be selected from \"dicot\", \"beta13\" ")
 
  if(length(diagNZ.strength) > 1 )  diagNZ.strength <- diagNZ.strength[1]
  if(diagNZ.strength< 0 | diagNZ.strength > 0.8)  stop("diagNZ.strength must be between zero and 0.8")
  
  ## exactZeroTh, mixProb, perturb.clust, mu, prob.sign, seed
  if(length(mixProb) > 1)  mixProb <- mixProb[1]
  if(mixProb < 0 | mixProb > 1)  stop("mixProb must be between zero and one")

  if(length(exactZeroTh) > 1)  exactZeroTh <- exactZeroTh[1]
  if(exactZeroTh < 0 )  stop("exactZeroTh must be larger or equal to zero")

  if(length(perturb.clust) > 1)  perturb.clust <- perturb.clust[1]
  if(perturb.clust < 0 | perturb.clust > 1)  stop("perturb.clust must be between zero and one")

  if(length(pdiff) > 1)  pdiff <- pdiff[1]
  if(pdiff < 0 | pdiff > 1)  stop("pdiff must be between zero and one")
  
  if(length(mu) == 1) mu <- rep(mu, sum(nnodesxcluster))
  if(length(mu) != sum(nnodesxcluster)) stop("length of mu must coincide with the sum of nnodesxcluster")
 
  if(length(probSign) > 1)  probSign <- probSign[1]
  if(probSign < 0 | probSign > 1)  stop("probSign must be between zero and one")
  
  if( diffType == "random" & length(seed) < nclusters) stop("seed must be a vector of size nclusters")
  if( (diffType == "cluster" | diffType == "mixed") & length(seed) < nclusters + 2) stop("seed must be a vector of size nclusters + 2")
  
   ## degree.hubs
   if(length(degree.hubs)==1)  degree.hubs <- rep(degree.hubs, nclusters)  
   if(any(degree.hubs <= 1))  stop("degree.hubs must be larger than one")

    if(length(nOtherEdges)==1)  nOtherEdges <- rep(nOtherEdges,nclusters)  
    if(any(nOtherEdges < 0))  stop("nOtherEdges must be larger than zero")

  listAll <- list(nobs = nobs, nclusters = nclusters, nnodesxcluster = nnodesxcluster,
  			      pattern = pattern, diffType = diffType, dataDepend = dataDepend, 
  			      low.strength = low.strength, sup.strength = sup.strength, 
  			      nhubs = nhubs, degree.hubs = degree.hubs, nOtherEdges = nOtherEdges, 
                  alpha = alpha, plus = plus, pdiff = pdiff, prob = prob, diagCCtype = diagCCtype,
                  diagNZ.strength = diagNZ.strength, mixProb = mixProb, perturb.clust = perturb.clust, 
                  mu = mu, probSign = probSign, seed = seed, exactZeroTh = exactZeroTh) 
  return(listAll)
}
