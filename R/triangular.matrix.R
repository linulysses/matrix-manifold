#'
#' log-cholesky metric on the Cholesky factor space
#' @param p the base point
#' @param u a tangent vector at p
#' @param v a tangent vector at p
#' @keywords internal
#' @export
rie.metric.ltpd_LogCholesky <- function(mfd,p,u,v,...)
{
    return( frobenius(lower.part(u),lower.part(v)) + sum(diag(u)*diag(v)/(diag(p))^2) )
}

#' geodesic on the Cholesky factor space
#' @param p starting point
#' @param u direction
#' @param t time
#' @keywords internal
#' @export
geodesic.ltpd_LogCholesky <- function(mfd,p,u,t,...)
{
    if(length(t) > 1)
    {
        res <- array(0,dim=c(dim(u),length(t)))
        for(i in 1:length(t))
        {
            s <- t[i]
            R <- lower.part(p) + s*lower.part(u)
            diag(R) <- diag(p) * exp(s*diag(u)/diag(p))
            res[,,i] <- R
        }
        
        return(res)
    }
    else
    {
        R <- lower.part(p) + t*lower.part(u)
        diag(R) <- diag(p) * exp(t*diag(u)/diag(p))
        return(R)
    }
    
}

#' log-cholesky distance on the Cholesky factor space
#' @param p a LTPD
#' @param q a LTPD
#' @keywords internal
#' @export
geo.dist.ltpd_LogCholesky <- function(mfd,p,q,...)
{
    if(is.matrix(p))
    {
        R <- frobenius(lower.part(p)-lower.part(q))
        R <- R + sum((log(diag(p))-log(diag(q)))^2)
        return(sqrt(R))
    }
    else # vector form
    {
        n <- length(p)
        d <- round(-0.5+sqrt(2*length(p)+0.25))
        R <- sum(((p[(d+1):n])-(q[(d+1):n]))^2)
        R <- R + sum((log((p[1:d]))-log((q[1:d])))^2)
        return(sqrt(R))
    }
    
}

#' log-cholesky log map on the Cholesky factor space
#' @param p the foot point of the log map
#' @param q the point to be mapped
#' @keywords internal
#' @export
rie.log.ltpd_LogCholesky <- function(mfd,p,q,...)
{
    R <- lower.part(q) - lower.part(p)
    diag(R) <- diag(p) * (log(diag(q))-log(diag(p)))
    return(R)
}

#' log-cholesky exponential map on the Cholesky factor space
#' @param p the base point
#' @param v a tangent vector at p
#' @keywords internal
#' @export
rie.exp.ltpd_LogCholesky <- function(mfd,p,v,...)
{
    return(geodesic.chol(p,v,1))
}



#' Frechet mean of positive lower-triangular matrix under Log-Cholesky metric
#' @param mfd an object created by \code{matrix.manifold}
#' @param S an \code{M*N*n} array of matrices or a list of \code{n} matrices, where \code{n} is the number of matrices
#' @return the Frechet mean of the matrices in \code{S}
#' @keywords internal
#' @export
frechet.mean.ltpd_LogCholesky <- function(mfd,S)
{
    R <- 0
    D <- 0
    n <- dim(S)[3]
    for(i in 1:n)
    {
        C <- as.matrix(S[,,i])
        R <- R + lower.part(C,T)
        D <- D + log(diag(C))
    }
    R <- R / n
    D <- D / n
    L <- R + diag(exp(D),nrow=mfd$dim[1],ncol=mfd$dim[2])
    return(L)
}

#' Parallel transport of a tangent vector from a point to another
#' @param p a LTPD matrix
#' @param q a LTPD matrix
#' @param v a tangent vector at \code{p}
#' @keywords internal
#' @export
parallel.transport.ltpd_LogCholesky <- function(mfd,p,q,v,...)
{
    Y <- lower.part(v,T) + diag(as.vector(diag(q)*diag(p)^(-1)*diag(v)),nrow=mfd$dim[1],ncol=mfd$dim[2])
    return(Y)
}


#' Generate a set of normally distributed random matrices that represent tangent vectors at some point. 
#' @param mfd an object created by \code{matrix.manifold}
#' @param n sample size
#' @param sig the standard deviation of the normal distribution
#' @param drop whether return the result as a matrix when \code{n=1}
#' @return an \code{M*N*n} array of \code{n} matrices, where \code{M*N} is the dimensions of matrices
#' @keywords internal
#' @export
rtvecor.ltpd_LogCholesky <- function(mfd,n=1,sig=1,drop=T)
{
    if(mfd$dim[1] > 1) mu <- diag(rep(0,mfd$dim[1]))
    else mu <- matrix(0,1,1)
    rltmatrix(mfd$dim[1],n=n,strict=F,mu=mu,
              sig=sig,dist='Gaussian',drop=drop)
}

#' Generate a set of random matrices on a matrix manifold
#' @param mfd an object created by \code{matrix.manifold}
#' @param n sample size
#' @param mu the Frechet mean. If \code{NULL} is given, then it is the identity element, e.g., for SPD, it is the identity matrix
#' @param sig the standard deviation of the normal distribution
#' @param drop whether return the result as a matrix when \code{n=1}
#' @return an \code{M*N*n} array of \code{n} matrices, where \code{M*N} is the dimensions of matrices
#' @details The generated samples have Frechet mean \code{mu}. The logarithmic maps of these samples at \code{mu} follow a isotropic D-dimensional normal distribution with isotropic variance \code{sig}, where D is the intrinsic dimension of the matrix manifold
#' @keywords internal
#' @export
rmatrix.ltpd_LogCholesky <- function(mfd,n=1,mu=NULL,sig=1,drop=T)
{
    
    if(is.null(mu)) mu <- diag(rep(1,mfd$dim[1]))
    stopifnot(is.ltpd(mu))
    
    S <- rtvecor(mfd,n,sig,drop=F)
    
    R <- array(0,c(mfd$dim[1],mfd$dim[1],n))
    for(i in 1:n)
    {
        R[,,i] <- rie.exp(mfd,mu,as.matrix(S[,,i]))
    }
    
    if(n==1 && drop){
        if(mfd$dim[1] > 1) return(R[,,1])
        else return(matrix(R[,,1],1,1))
    }
    else return(R)
}



#' generate a set of random lower triangular matrix of dimension d
#' @keywords internal
#' @importFrom stats runif rnorm
rltmatrix <- function(d,n=1,strict=F,mu=NULL,sig=1,dist='Gaussian',drop=T)
{
    if(is.null(mu)) mu <- diag(d)
    
    if(d==1 && strict==T) stop('strict lower triangular matrix does exist when d=1')
    S <- array(0,c(d,d,n))
    for(i in 1:n)
    {
        if(dist == 'uniform')
            S[,,i] <- lower.part(matrix(runif(d*d),d,d),strict)
        else if(dist == 'Gaussian')
        {
            if(d==1)
            {
                S[,,i] <- exp(rnorm(1,mean=mu,sd=sig))
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
    
    if(n==1 && drop) 
    {
        if(d > 1) return(S[,,1])
        else return(matrix(S[,,1],1,1))
    }
        
    else return(S)
}



#' generate a random upper triangular matrix of dimension d
#' @keywords internal
rutmatrix <- function(d,n=1,strict=F,mu=diag(d),sig=1,dist='Gaussian',drop=T)
{
    batch.t(rltmatrix(d,n,strict,t(mu),sig,dist,drop))
}

