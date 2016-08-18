softA <- function (a, lam, penalize.diagonal) 
{
    out <- sign(a) * pmax(0, as.matrix(abs(a) - lam))
    if (!penalize.diagonal) 
        diag(out) <- diag(as.matrix(a))
    return(out)
}


