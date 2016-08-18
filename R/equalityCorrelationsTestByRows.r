eqCorTestByRows <- function(D1, D2 = NULL, testStatistic = c("AS", "max"), nite = 200, 
							paired = FALSE, exact = TRUE, subMatComp = FALSE, iniP = 1, 
							finP = NULL, conf.level = 0.95)
{ 
  
  ## checks
  testStatistic <- controlseqCorTestByRows(D1 = D1, D2 = D2, testStatistic = testStatistic, 
  		nite = nite, paired = paired, conf.level = conf.level)

  if(!is.null(D2))
  {
   if(is.null(finP)) finP <- dim(D1)[2]
   
   if(iniP - finP > -1 & subMatComp)
    stop("iniP must be smaller than finP ")

   if(finP > dim(D1)[2] & subMatComp)
    stop("finP must be smaller than the dimension")
    
   if(iniP <=0 & subMatComp)
    stop("iniP must be larger than zero")
  }
  
  ## Functions
  if(is.null(D2))
   obj <- correlationTest.byRows(D1, testStatistic = testStatistic, nite = nite,
   								 subMatComp = subMatComp, conf.level = conf.level)
  else
   obj <- eqCorrelationTest.byRows(D1 = D1, D2 = D2, testStatistic = testStatistic,
  								   nite = nite, paired = paired, exact = exact,
  								   subMatComp = subMatComp, iniP = iniP, finP = finP,
  								   conf.level = conf.level)
  
  class(obj) <- "eqCorTestByRows"
  return(obj)
}
