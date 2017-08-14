 
########################## Partial correlation and data simulator for paired data ##############
pcorSimulatorTwo <- function(nobs, nclusters, nnodesxcluster, pattern="hubs", diffType=c("random","cluster", "mixed"),
                               nhubs=5, degree.hubs = 20, low.strength = 0.5, sup.strength = 0.9, 
                               nOtherEdges = 30, alpha = 2.3, plus=0, seed=sample(1:10000,nclusters+2), 
                               prob = 0.05, perturb.clust=0, mu = rep(0, sum(nnodesxcluster)), pdiff=0,  mixProb=0.5, probSign = 0.5,  
                               exactZeroTh = 0.05)
{
      
   ## Initial partial correlation matrix for random differential matrices
   if(pdiff<=0||diffType=="random"){
       PCOR <-  pcorSimulator(nobs = nobs, nclusters = nclusters, nnodesxcluster = nnodesxcluster, 
                       pattern = pattern, nhubs = nhubs, degree.hubs = degree.hubs, 
                       low.strength = low.strength, sup.strength=sup.strength, nOtherEdges = nOtherEdges, 
                       alpha = alpha, plus=plus, seed = seed, prob=prob, perturb.clust = perturb.clust, 
                       mu = mu, probSign=probSign)
                           
       PREC.MAT.0 	<- PCOR$omega
       PREC.MAT.0 	<- cov2cor(PREC.MAT.0)
       PREC.MAT.1 	<- PREC.MAT.0 
       P 			<- dim(PREC.MAT.1)[1]
   }
     
   ## Initial partial correlation matrix for cluster or mixed differential matrices   
   if(pdiff > 0 &&(diffType == "cluster" | diffType=="mixed")){
       dnod    <- round((pdiff*sum(nnodesxcluster)/2))
       if(diffType=="mixed") dnod <- round(dnod*(mixProb))

       if(dnod > 4){
         if(pattern=="hubs" & diffType=="cluster" &length(degree.hubs)== nclusters) 
           degree.hubs <- c(degree.hubs,round(rep(degree.hubs[1] *(dnod/nnodesxcluster[1]),2)))
         if(pattern=="hubs" & diffType=="cluster" &length(nOtherEdges)== nclusters) 
           nOtherEdges <- c(nOtherEdges,round(rep(nOtherEdges[1] *(2^dnod/2^nnodesxcluster[1]),2)))
           
         PCOR       <- pcorSimulator(nobs = nobs, nclusters = nclusters + 2, nnodesxcluster = c(nnodesxcluster,dnod,dnod), 
                       pattern = pattern, nhubs = nhubs, degree.hubs = degree.hubs, 
                       low.strength=low.strength, sup.strength=sup.strength, nOtherEdges = nOtherEdges, 
                       alpha = alpha, plus=plus, seed = seed, prob=prob, perturb.clust = perturb.clust, 
                       mu = c(mu, rep(0,dnod*2)), probSign=probSign)
       
        PREC.MAT.0 	<- PCOR$omega
        PREC.MAT.0 	<- cov2cor(PREC.MAT.0)
        PREC.MAT.1 	<- PREC.MAT.0 
        P 			<- dim(PREC.MAT.1)[1]
       
        PREC.MAT.0[(sum(nnodesxcluster)+1):(sum(nnodesxcluster)+dnod),(sum(nnodesxcluster)+1):(sum(nnodesxcluster)+dnod)]<- 0
        PREC.MAT.1[(sum(nnodesxcluster)+dnod+1):(sum(nnodesxcluster)+dnod*2),(sum(nnodesxcluster)+dnod+1):(sum(nnodesxcluster)+dnod*2)]<- 0
        diag(PREC.MAT.0) <- 1
        diag(PREC.MAT.1) <- 1
        PREC.MAT.0       <- as.matrix(PREC.MAT.0)
        PREC.MAT.1       <- as.matrix(PREC.MAT.1)
        nclusters        <- nclusters+2
        nnodesxcluster   <- c(nnodesxcluster,dnod,dnod)
       }
       else{
        PCOR       <- pcorSimulator(nobs = nobs, nclusters = nclusters, nnodesxcluster = nnodesxcluster, 
                       pattern = pattern, nhubs = nhubs, degree.hubs = degree.hubs, 
                       low.strength=low.strength, sup.strength=sup.strength, nOtherEdges = nOtherEdges, 
                       alpha = alpha, plus=plus, seed = seed, prob=prob, perturb.clust = perturb.clust, 
                       mu = mu, probSign=probSign)
       
        PREC.MAT.0 	<- PCOR$omega
        PREC.MAT.0 	<- cov2cor(PREC.MAT.0)
        PREC.MAT.1 	<- PREC.MAT.0 
        P 			<- dim(PREC.MAT.1)[1]
       }
   }
 
  ## Differential partial correlation matrix for random and mixed   
  if(pdiff > 0 && (diffType == "random" | diffType == "mixed")){
    CORR        <- cov2cor(solve(PREC.MAT.0))
    PREC.MAT.0  <- pseudoinverse(CORR)
    set.seed(seed[1]*3)
    id <- runif(P,0,1) < pdiff
    CORR[id, !id] <- 0
    CORR[!id, id] <- 0
    PREC.MAT.1  <- pseudoinverse(CORR)
   
    PREC.MAT.0[PREC.MAT.0<exactZeroTh] <- 0
    PREC.MAT.1[PREC.MAT.1<exactZeroTh] <- 0
  }
     
  return(list(omega0 = PREC.MAT.0, omega1 = PREC.MAT.1, nclusters=nclusters, nnodesxcluster=nnodesxcluster ))

}


