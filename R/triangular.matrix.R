#' transform a (set of) lower triangular matrix into a vector representation, diagonal by diagonal, starting from the main diagonal
#' @keywords internal
lower.to.vec.atomic <- function(L)
{
    if(is.vec(L) && length(L)==1) return(L)
    
    if(!is.matrix(L)) stop('L must be a matrix object')
    
    d <- dim(L)[1]
    v <- rep(0,d*(d+1)/2)
    v[1:d] <- diag(L)
    p <- d
    for(i in 2:d)
        for(j in 1:(d-i+1))
        {
            p <- p + 1
            v[p] <- L[i+j-1,j]
        }
    
    attr(v,'format') <- 'lower.in.vec'
    return(v)        
}

#' transform a (set of) lower triangular matrix into a vector representation, diagonal by diagonal, starting from the main diagonal
#' @keywords internal
lower.to.vec <- function(L,drop=F)
{
    if(is.vec(L) && length(L)==1) return(L)
    if(is.matrix(L)) L <- array(L,c(dim(L),1))
    
    if(is.array(L))
    {
        dL <- dim(L)
        X <- matrix(0,dL[2]*(dL[2]+1)/2,dL[3])
        for(i in 1:dL[3])
        {
            X[,i] <- lower.to.vec.atomic(L[,,i])
        }
        if(drop && dim(L)[3]==1) X <- X[,1]
        
        attr(X,'format') <- 'lower.in.vec'
        return(X)
    }
    else if(is.list(L))
    {
        lapply(L,lower.to.vec.atomic)
    }
}

#' recover a lower triangular matrix from its vector representation
#' @keywords internal
vec.to.lower.atomic <- function(v)
{
    d <- round(-0.5+sqrt(2*length(v)+0.25))
    L <- matrix(0,d,d)
    p <- 0
    for(i in 1:d)
        for(j in 1:(d-i+1))
        {
            p <- p + 1
            L[i+j-1,j] <- v[p]
        }
    attr(L,'format') <- 'lower.in.matrix'
    return(L)         
}

#' recover a (set of) lower triangular matrix from its vector representation
#' @keywords internal
vec.to.lower <- function(v,drop=F)
{
    if(is.vec(v))
    {
        v <- as.matrix(v)
        attr(v,'format') <- 'lower.in.vec'
    } 
    
    if(is.matrix(v))
    {
        stopifnot(attr(v,'format') == 'lower.in.vec')
        
        n <- ncol(v)
        d <- round(-0.5+sqrt(2*length(v[,1])+0.25))
        L <- array(0,c(d,d,n))
        for(i in 1:n)
        {
            L[,,i] <- vec.to.lower.atomic(v[,i])
        }
        if(drop && n == 1) return(L[,,1])
        else return(L)
    }
    else if(is.list(v))
    {
        lapply(v,function(u) vec.to.lower.atomic(v[[i]]))
    }
}
