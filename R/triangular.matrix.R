


#' generate a set of random lower triangular matrix of dimension d
#' @keywords internal
gen.lower <- function(d,n=1,strict=F,mu=diag(d),sig=1,dist='Gaussian',drop=T)
{
    S <- array(0,c(d,d,n))
    for(i in 1:n)
    {
        if(dist == 'uniform')
            S[,,i] <- lower.part(matrix(runif(d*d),d,d),strict)
        else if(dist == 'Gaussian')
        {
            if(d==1)
            {
                S[,,i] <- exp(rnorm(1,mean=log(mu),sd=sig))
            }
            else
            {
                v.mu <- lower.to.vec(mu)[-(1:d)]
                Sig <- sig * diag(d*(d-1)/2)
                L <- vec.to.lower(c(rep(0,d),MASS::mvrnorm(1,mu=v.mu,Sigma=Sig)),drop=T)
                if(!strict)
                {
                    diag(L) <- exp(MASS::mvrnorm(1,mu=log(diag(mu)),Sigma=sig*diag(d)))
                }
                S[,,i] <- L
            }
        }
    }
    
    if(n==1 && drop) return(S[,,i])
    else return(S)
}



#' generate a random upper triangular matrix of dimension d
#' @keywords internal
gen.upper <- function(d,n=1,strict=F,mu=diag(d),sig=1,dist='Gaussian',drop=T)
{
    batch.t(gen.lower(d,n,strict,t(mu),sig,dist,drop))
}


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
