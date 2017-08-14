
lambdaSelection <- function(obj, criterion=c("PC", "AGNES", "A-MSE", "VUL", "STARS", "AIC", "BIC", "eBIC"), ...)
{
  if(class(obj) != "huge" & class(obj) != "tiger" & class(obj) != "wfgl")
    stop("object has to be of class huge, tiger or wfgl")
  if(class(obj) == "wfgl")
  	obj$lambda <- obj$lambda1
  
  if(length(obj$lambda) < 5)
    stop("lambda sequence has to be of length 5 or higher")
  
  pcriterion <- c("PC","AGNES","A-MSE","STARS","VUL", "BIC", "AIC", "eBIC")
   if(length(criterion ) > 1)
    warning("criterion attribute length is larger than one. Only the first component will be used")
 
  criterion <- pcriterion[pmatch(criterion[1], pcriterion)]
  if (is.na(criterion)) stop("criterion is not well define. It must be selected from \"PC\", \"AGNES\", \"A-MSE\", \"VUL\" or \"STARS\" ")
  
  
  if(criterion[1] == "PC")
    ret.list <- pcLambdaSelection(obj)
  if(criterion[1] == "AGNES")
    ret.list <- agnesLambdaSelection(obj, ...)
  if(criterion[1] == "A-MSE")
    ret.list <- amseLambdaSelection(obj, ...)
  if(criterion[1] == "STARS")
    ret.list <- huge.select(obj, criterion = "stars", ...)
  if(criterion[1] == "VUL")
    ret.list <- vulLambdaSelection(obj, ...)
  if(criterion[1] == "AIC"|criterion[1] == "BIC"|criterion[1] == "eBIC")
    ret.list <- aicAndbicLambdaSelection(obj, criterion = criterion[1], ...)
  return(ret.list)
}