controlseqCorTestByRows <- function(D1, D2, testStatistic, nite, paired, conf.level)
{
  ## dimensions
  if(!is.null(D2)){
   if(dim(D1)[2] != dim(D2)[2])
    stop("Datasets dimensions must be equal")
  
   if(dim(D1)[1] != dim(D2)[1] & paired)
    stop("Datasets sample sizes must be equal when paired = TRUE")
  } 
  ## prospective testStatistic
  ptestStatistic <- c("AS", "max", "exc")
  testStatistic2 <- match.arg(testStatistic, ptestStatistic, several.ok = TRUE)
  if(length(testStatistic2) <  length(testStatistic) ){
   if(length(testStatistic2) >= 1)
    warning("not all elements in testStatistic are well defined")
  }
  testStatistic <- testStatistic2   
 
  ## nite, iniP, finP, conf.level
  if(nite < 10)
    stop("not enought permutations. Increase nite to be larger than 10")
    
  if(conf.level < 0.5 | conf.level >= 1)
    stop("conf.level must be between [0.5,1)")
  
  return(testStatistic)
}
