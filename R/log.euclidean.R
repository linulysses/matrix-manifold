#' Differential of matrix log of SPD manifold
#' @param S base point
#' @param V tangent vector
#' @keywords internal
d.logm <- function(S,V)
{
    eig <- eigen(S,symmetric = T)
    eigvals <- eig$values
    eigvecs <- eig$vectors
    
    d <- dim(S)[1]
    
    log.eigvals <- log(eigvals)
    denom <- matrix(eigvals,d,d,byrow=F)
    denom <- denom - t(denom)
    numer <- matrix(log.eigvals,d,d,byrow=F)
    numer <- numer - t(numer)
    
    numer <- ifelse(denom==0,yes=matrix(1,d,d),no=numer)
    denom <- ifelse(denom==0,yes=matrix(eigvals,d,d,byrow=F),no=denom)
    
    tmp <- t(eigvecs) %*% V %*% eigvecs
    
    R <- eigvecs %*% ((numer/denom) * tmp) %*% t(eigvecs)
    
    return(R)
}

#' Differential of matrix exponential of SPD manifold
#' @param S base point
#' @param V tangent vector
#' @keywords internal
d.expm <- function(S,V)
{
    eig <- eigen(S,symmetric = T)
    eigvals <- eig$values
    eigvecs <- eig$vectors
    
    d <- dim(S)[1]
    
    exp.eigvals <- exp(eigvals)
    denom <- matrix(eigvals,d,d,byrow=F)
    denom <- denom - t(denom)
    numer <- matrix(exp.eigvals,d,d,byrow=F)
    numer <- numer - t(numer)
    
    numer <- ifelse(denom==0,yes=matrix(exp.eigvals,d,d,byrow=F),no=numer)
    denom <- ifelse(denom==0,yes=matrix(1,d,d),no=denom)
    
    tmp <- t(eigvecs) %*% V %*% eigvecs
    
    R <- eigvecs %*% ((numer/denom) * tmp) %*% t(eigvecs)
    
    return(R)
}

#' Log-Euclidean group multiplication
#' @keywords internal
mul.LE <- function(P,Q)
{
    expm::expm(expm::logm(P)+expm::logm(Q))
}


#' Compute the Riemannian exponential map
#' @param mfd the manifold object created by \code{matrix.manifold}
#' @param p the foot point of the exponential map
#' @param v a tangent vector at \code{p}
#' @param ... other parameters (primarily for computation speedup)
#' @return a matrix that is the exponential map of \code{v}
#' @export
rie.exp.spd_LogEuclidean <- function(mfd,p,v,...)
{
    make.sym(expm::expm(expm::logm(p)+d.logm(p,v)))
}

#' Compute the Riemannian logarithm map
#' @param mfd the manifold object created by \code{matrix.manifold}
#' @param p the foot point of the logarithm map
#' @param q a point on the manifold
#' @param ... other parameters (primarily for computation speedup)
#' @return a matrix that is the logarithm map of \code{q}
#' @export
rie.log.spd_LogEuclidean <- function(mfd,p,q,...)
{
    logp <- expm::logm(p)
    make.sym(d.expm(logp, expm::logm(q) - logp))
}

#' Compute the Riemannian metric 
#' @param mfd the manifold object created by \code{matrix.manifold}
#' @param p the foot point of the metric
#' @param u a tangent vector at \code{p}
#' @param v a tangent vector at \code{q}
#' @param ... other parameters (primarily for computation speedup)
#' @return a number that is the Riemannian metric of \code{u} and \code{v}
#' @export
rie.metric.spd_LogEuclidean <- function(mfd,p,u,v,...)
{
    sum(d.logm(p,u)*d.logm(p,v))
}

#' Compute the geodesic starting at some point on a manifold
#' @param mfd the manifold object created by \code{matrix.manifold}
#' @param p the starting point of the geodesic
#' @param u a tangent vector at \code{p}, representing the direction of the geodesic
#' @param t a vector of nonnegative real numbers
#' @param ... other parameters (primarily for computation speedup)
#' @return an array of matrices which are the geodesic evaluated at \code{t}. The last dimension of the array corresponds to \code{t}
#' @export
geodesic.spd_LogEuclidean <- function(mfd,p,u,t,...)
{
    d <- mfd$dim[1]
    if(length(t) > 1)
    {
        R <- array(0,c(d,d,length(t)))
        for(i in 1:length(t))
            R[,,i] <- make.sym(rie.exp(mfd,p,t[i]*u))
    }
    else
    {
        R <- make.sym(rie.exp(mfd,p,t*u))
    }
    
    return(R)
}

#' Compute the geodesic distance
#' @param mfd the manifold object created by \code{matrix.manifold}
#' @param p a point on the manifold
#' @param q a point on the manifold
#' @param ... other parameters (primarily for computation speedup)
#' @return the geodesic distance between \code{p} and \code{q}
#' @export
geo.dist.spd_LogEuclidean <- function(mfd,p,q,...)
{
    sqrt(sum((expm::logm(p)-expm::logm(q))^2))
}

#' Parallel transport of a tangent vector from one point to another along the geodesic between the two points
#' @param mfd the manifold object created by \code{matrix.manifold}
#' @param p a point on the manifold
#' @param q a point on the manifold
#' @param v a tangent vector at \code{p}
#' @param ... other parameters (primarily for computation speedup)
#' @return a tangent vector at \code{q}
#' @export
parallel.transport.spd_LogEuclidean <- function(mfd,p,q,v,...) 
{
    S <- mul.LE(q,inv(p))
    d.expm(expm::logm(S)+expm::logm(p), d.logm(p,v))
}

#' Generate a set of random matrices that represent tangent vectors at some point. 
#' @param mfd an object created by \code{matrix.manifold}
#' @param n sample size
#' @param sig the standard deviation of the normal distribution
#' @param drop whether return the result as a matrix when \code{n=1}
#' @return an \code{M*N*n} array of \code{n} matrices, where \code{M*N} is the dimensions of matrices
#' @details The generated samples have expectation zero and follow a isotropic D-dimensional normal distribution with isotropic variance \code{sig}, where D is the intrinsic dimension of the matrix manifold
#' @export
rtvecor.spd_LogEuclidean <- function(mfd,n=1,sig=1,drop=T)
{
    return(rsym(d=mfd$dim[1],n=n,sig=sig,drop=drop))
}

#' Generate a set of random matrices on a matrix manifold
#' @param mfd an object created by \code{matrix.manifold}
#' @param n sample size
#' @param mu the Frechet mean. If \code{NULL} is given, then it is the identity element, e.g., for SPD, it is the identity matrix
#' @param sig the standard deviation of the normal distribution
#' @param drop whether return the result as a matrix when \code{n=1}
#' @return an \code{M*N*n} array of \code{n} matrices, where \code{M*N} is the dimensions of matrices
#' @details The generated samples have Frechet mean \code{mu}. The logarithmic maps of these samples at \code{mu} follow a isotropic D-dimensional normal distribution with isotropic variance \code{sig}, where D is the intrinsic dimension of the matrix manifold
#' @export
rmatrix.spd_LogEuclidean <- function(mfd,n=1,mu=NULL,sig=1,drop=T)
{
    if(is.null(mu)) mu <- diag(rep(1,mfd$dim[1]))
    stopifnot(is.spd(mu))
    
    S <- rtvecor(mfd,n,sig,drop=F)
    
    R <- array(0,c(mfd$dim[1],mfd$dim[1],n))
    for(i in 1:n)
    {
        R[,,i] <- rie.exp(mfd,mu,S[,,i])
    }
    
    if(n==1 && drop) return(as.matrix(R[,,1]))
    else return(R)
}

#' Frechet mean of matrices
#' @param mfd an object created by \code{matrix.manifold}
#' @param S an \code{M*N*n} array of matrices or a list of \code{n} matrices, where \code{n} is the number of matrices
#' @return the Frechet mean of the matrices in \code{S}
#' @export
frechet.mean.spd_LogEuclidean <- function(mfd,S)
{
    R <- 0
    if(is.list(S))
    {
        for(i in 1:length(S))
        {
            R <- R + expm::logm(S[[i]])
        }
        R <- R / length(S)
        return(expm::expm(R))
    }
    else if(is.array(S))
    {
        n <- dim(S)[3]
        d <- dim(S)[1]
        for(i in 1:n)
        {
            R <- R + expm::logm(as.matrix(S[,,i]))
        }
        R <- R/n
        return(expm::expm(R))
    }
    else if(is.matrix(S)) return(S)
    else stop('S must be an array, a list or a matrix')
    
}