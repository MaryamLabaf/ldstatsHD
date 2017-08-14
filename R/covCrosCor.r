covCrosCor <- function(rhojk, rhohm, rhokh, rhokm, rhojh, rhojm)
{
        rhokj <- rhojk
        rhomh <- rhohm
        AU <- 1/(2*(1-rhojk^2)*(1-rhohm^2)) * (
          (rhojh - rhojk*rhokh)*(rhokm - rhokh*rhohm) +
          (rhojm - rhojh*rhohm)*(rhokh - rhokj*rhojh) +
          (rhojh - rhojm*rhomh)*(rhokm - rhokj*rhojm) +
          (rhojm - rhojk*rhokm)*(rhokh - rhokm*rhomh))
        return(AU)
}

