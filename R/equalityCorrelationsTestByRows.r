eqCorTestByRows <- function(D1, D2 = NULL, testStatistic = c("AS", "max"), nite = 200, 
							paired = FALSE, exact = TRUE, whichRows = NULL, conf.level = 0.95,...)
{ 
  
  ## checks
  testStatistic <- controlseqCorTestByRows(D1 = D1, D2 = D2, testStatistic = testStatistic, 
  		nite = nite, paired = paired, conf.level = conf.level)
  		  
  ## Functions
  if(is.null(D2))
   obj <- correlationTest.byRows(D1, testStatistic = testStatistic, nite = nite,
   								 whichRows = whichRows, conf.level = conf.level,...)
  else
   obj <- eqCorrelationTest.byRows(D1 = D1, D2 = D2, testStatistic = testStatistic,
  								   nite = nite, paired = paired, exact = exact,
  								   whichRows = whichRows,conf.level = conf.level,...)
  
  class(obj) <- "eqCorTestByRows"
  return(obj)
}
