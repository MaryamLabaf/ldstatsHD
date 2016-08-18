estradaIndex <-function(A)
{
      EIG                   <- eigen(A)
      return( sum(exp(EIG$val)))
}
