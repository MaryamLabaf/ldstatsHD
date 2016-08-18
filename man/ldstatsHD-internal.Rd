
\name{ldstatsHD-inl}
\alias{ldstatsHD-internal}

\alias{covCrosCor}
\alias{ztransf}
\alias{rowVars}
\alias{rowScale}
\alias{softA}
\alias{lowerTri}
\alias{upperTri}
\alias{lowerTriInd}
\alias{lowerTriMatInd}
\alias{pvaluesTriangleWeakestEdge}
\alias{asMatrix2}
\alias{findTriangles}
\alias{findAlmTriangles}
\alias{graphStructure}
\alias{connBtwClusters}
\alias{controlsPcorSimulator}
\alias{controlsPcorSimulatorJoint}
\alias{ztransfCorrDiff}
\alias{thetaFromCor}
\alias{ztransfCorrDiffsub}
\alias{meanCors}
\alias{controlseqCorTestByRows}
\alias{controlwfgl}
\alias{controlwfrl}
\alias{flsaFast}
\alias{pcor2jointGraph}
\alias{squareWeight}
\alias{wfglaux}
\alias{wfrlaux}
\alias{findingLambda2}

\title{Internal  functions}
\description{Internal  functions}
\usage{
covCrosCor(rhojk, rhohm, rhokh, rhokm, rhojh, rhojm)
ztransf(r)
rowVars(mat)
rowScale(mat)
softA(a, lam, penalize.diagonal) 
lowerTri(X)
upperTri(X)
lowerTriInd(p)
lowerTriMatInd(x, P)
pvaluesTriangleWeakestEdge(A1, A2, alpha=0.01, D1, D2, equalBetweenCond = TRUE)
asMatrix2(M)
findTriangles(A1)
findAlmTriangles(A1)
graphStructure(nnodes, low.strength, sup.strength, pattern ="hubs", nhubs, 
               alpha = 2.3, degree.hubs, nOtherEdges, orderNodes=0, 
               seed = 200, plus = 0, prob=0.05, probSign=0.5)
connBtwClusters(PREC.MAT, perturb.clust = 0)
controlsPcorSimulator(nobs, nclusters, nnodesxcluster, pattern, 
                      low.strength, sup.strength, nhubs, degree.hubs, nOtherEdges, 
                      alpha, plus, prob, perturb.clust, mu, probSign, seed)
controlsPcorSimulatorJoint(nobs, nclusters, nnodesxcluster, pattern, 
                      diffType, dataDepend, low.strength, sup.strength, pdiff, nhubs, 
                      degree.hubs, nOtherEdges, alpha, plus, prob, perturb.clust, mu, 
                      diagCCtype, diagNZ.strength, mixProb, probSign, seed, exactZeroTh)
ztransfCorrDiff(D1, D2, thetaKnown=NULL, biasCorr=NULL )
thetaFromCor(COR, P)
ztransfCorrDiffsub(D1, D2, row, thetaKnown = NULL,  diagCor )
meanCors(mat, scaleT = TRUE, sampleFrom = FALSE, nitSamp = 1000)
controlseqCorTestByRows(D1, D2, testStatistic, nite, paired, conf.level)
controlwfgl(D1, D2, lambda1, lambda2, paired, automLambdas, 
	 		sigmaEstimate, pairedEst, maxiter, tol, nsubset, weights,
	 		rho, rho.increment, triangleCorrection, alphaTri, temporalFolders)
controlwfrl(D1, D2, lambda1, lambda2, paired, automLambdas, sigmaEstimate, 
            maxiter, tol, nsubset, rho, rho.increment, notOnlyLambda2)
flsaFast(A, L, lam1, lam2, corsOmega = NULL, P, ind12, 
        nsubset = 10000, automLambdas = TRUE, penalize.diagonal = FALSE, 
        sigmaEstimate = "CRmad", temporalFolders, notOnlyLambda2 = TRUE)
pcor2jointGraph(M1, M2, minn=0, col=c("blue","red","green"), vertex.size=2, 
  C1=NULL,C2=NULL, edgesThicknees=FALSE, zoomThick = 10, ...)
degrees(A)
squareWeight(at, th, weight1)
wfglaux(D1, D2, lambda1, lambda2, paired = TRUE, automLambdas = TRUE, 
sigmaEstimate = "CRmad", pairedEst = "Reg-based-sim", maxiter = 30, 
 tol = 1e-05, nsubset = 10000, weights = c(1,1), rho = 1, rho.increment = 1, 
 triangleCorrection = TRUE, alphaTri = 0.01, temporalFolders = FALSE, 
 theta0 = NULL, corEst2 = NULL, ind12 = NULL, S = NULL, notOnlyLambda2 = TRUE)	
wfrlaux(D1, D2, lambda1, lambda2, automLambdas=FALSE, paired = FALSE, 
  sigmaEstimate = "CRmad", maxiter = 30, tol=1e-05, nsubset = 10000,  
  rho = 1, rho.increment = 1, notOnlyLambda2 = TRUE, corsOmega = NULL, 
  thetae0 = NULL, theta = NULL, C12 = NULL, ind12 = NULL)
findingLambda2(n, lam1, lam2, tol, maxiter=100)



}

\details{These are not to be called by the user}
\keyword{internal}
