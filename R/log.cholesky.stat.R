
#' generate a set of random symmetric matrix of dimension d
#' @param dist the distribution, either 'uniform' or 'Gaussian'.
#' @param mu the mean
#' @param sig a scalar representing the dispersion around the mean
#' @param drop drop the last dimension if \code{n=1}
#' @keywords internal
gen.sym <- function(d,n=1,mu=matrix(0,d,d),sig=1,dist='Gaussian',drop=T)
{
    S <- array(0,c(d,d,n))
    for(i in 1:n)
    {
        if(dist == 'uniform')
        {
            S[,,i] <- make.sym(matrix(runif(d*d),d,d))
        }
        else if(dist == 'Gaussian')
        {
            if(d==1)
            {
                S[,,i] <- rnorm(1,mean=mu,sd=sig)
            }
            else
            {
                v.mu <- lower.to.vec(mu)[-(1:d)] * sqrt(2)
                Sig <- sig * diag(d*(d-1)/2)
                v <- c(rep(0,d),MASS::mvrnorm(1,mu=v.mu,Sigma=Sig))
                L <- vec.to.lower(v,drop=T) / sqrt(2)
                A <- L + t(L)
                diag(A) <- MASS::mvrnorm(1,mu=diag(mu),Sigma=sig*diag(d))
                S[,,i] <- A
            }
        }
    }
    
    if(n==1 && drop) S <- S[,,i]
    attr(S,'format') <- 'sym.in.matrix'
    return(S)
    
}


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

#' transpose all matrices in an array
#' @keywords internal
batch.t <- function(S)
{
    if(is.matrix(S)) return(t(S))
    else
    {
        d <- dim(S)
        R <- array(0,c(d[2],d[1],d[3]))
        for(i in 1:d[3])
        {
            R[,,i] <- t(S[,,i])
        }
    }
    return(R)
}


#' generate a random upper triangular matrix of dimension d
#' @keywords internal
gen.upper <- function(d,n=1,strict=F,mu=diag(d),sig=1,dist='Gaussian',drop=T)
{
    batch.t(gen.lower(d,n,strict,t(mu),sig,dist,drop))
}



#' Generate a set of normally distributed random matrices that represent tangent vectors at some point. 
#' @param mfd an object created by \code{create.matrix.manifold}
#' @param n sample size
#' @param sig the standard deviation of the normal distribution
#' @return an \code{M*N*n} array of \code{n} matrices, where \code{M*N} is the dimensions of matrices
#' @keywords internal
#' @export
gen.tangent.vectors.spd_LogCholesky <- function(mfd,n=1,sig=1,drop=T)
{
    return(gen.sym(d=mfd$dim[1],n=n,sig=sig,drop=drop))
}

#' Generate a set of random matrices on a matrix manifold
#' @param mfd an object created by \code{create.matrix.manifold}
#' @param n sample size
#' @param mu the Frechet mean. If \code{NULL} is given, then it is the identity element, e.g., for SPD, it is the identity matrix
#' @param sig the standard deviation of the normal distribution
#' @return an \code{M*N*n} array of \code{n} matrices, where \code{M*N} is the dimensions of matrices
#' @details The generated samples have Frechet mean \code{mu}. The logarithmic maps of these samples at \code{mu} follow a isotropic D-dimensional normal distribution with isotropic variance \code{sig}, where D is the intrinsic dimension of the matrix manifold
#' @keywords internal
#' @export
gen.matrices.spd_LogCholesky <- function(mfd,n=1,mu=NULL,sig=1,drop=T)
{

    if(is.null(mu)) mu <- diag(rep(1,mfd$dim[1]))
    stopifnot(is.spd(mu))
    
    S <- gen.tangent.vectors(mfd,n,sig,drop=F)
    
    R <- array(0,c(mfd$dim[1],mfd$dim[1],n))
    for(i in 1:n)
    {
        R[,,i] <- rie.exp(mfd,mu,S[,,i])
    }
    
    if(n==1 && drop) return(R[,,i])
    else return(R)
}

#' Frechet mean of matrices
#' @param mfd an object created by \code{create.matrix.manifold}
#' @param S an \code{M*N*n} array of matrices or a list of \code{n} matrices, where \code{n} is the number of matrices
#' @return the Frechet mean of the matrices in \code{S}
#' @keywords internal
#' @export
frechet.mean.spd_LogCholesky <- function(mfd,S)
{
    R <- 0
    D <- 0
    if(is.list(S))
    {
        for(i in 1:length(S))
        {
            C <- ch(S[[i]])
            R <- R + lower.part(C,T)
            D <- D + log(diag(C))
        }
        R <- R / length(S)
        D <- D / length(S)
        L <- R + diag(expm.diag(D))
        return(L %*% t(L))
    }
    else if(is.array(S))
    {
        n <- dim(S)[3]
        d <- dim(S)[1]
        for(i in 1:n)
        {
            C <- ch(S[,,i])
            R <- R + lower.part(C,T)
            D <- D + log(diag(C))
        }
        R <- R / n
        D <- D / n
        L <- R + diag(exp(D),nrow=d,ncol=d)
        return(L %*% t(L))
    }
    else if(is.matrix(S)) return(S)
    else stop('S must be an array, a list or a matrix')
    
}

#' Center a sample of matrices in Euclidean space
#' @param S an array of matrices where the last dimension corresponds to sample size
#' @param mu the center
#' @keywords internal
center.matrices <- function(S,mu=NULL)
{
    if(is.null(mu)) mu <- euclidean.mean(S)

    n <- dim(S)[3]
    for(i in 1:n)
    {
        S[,,i] <- S[,,i] - mu
        
    }
    return(S)
}

#' Euclidean mean of a sample of matrices
#' @param S an array of matrices where the last dimension corresponds to sample size
#' @keywords internal
euclidean.mean <- function(S)
{
    R <- 0
    n <- dim(S)[3]
    for(i in 1:n)
    {
        R <- R + S[,,i]
        
    }
    R <- R / n
    return(R)
}

#' Frechet mean of positive lower-triangular matrix under Log-Cholesky metric
#' @keywords internal
frechet.mean.chol <- function(S)
{
    R <- 0
    D <- 0
    n <- dim(S)[3]
    for(i in 1:n)
    {
        C <- S[,,i]
        R <- R + lower.part(C,T)
        D <- D + log(diag(C))
    }
    R <- R / n
    D <- D / n
    L <- R + diag(exp(D))
    return(L)
}
