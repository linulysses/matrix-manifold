#'
#' log-cholesky metric on the Cholesky factor space
#' @keywords internal
metric.chol <- function(L,X,Y)
{
    return( frobenius(lower.part(X),lower.part(Y)) + sum(diag(X)*diag(Y)/(diag(L))^2) )
}

#' geodesic on the Cholesky factor space
#' @keywords internal
geodesic.chol <- function(L,X,t)
{
    if(length(t) > 1)
    {
        res <- array(0,dim=c(dim(X),length(t)))
        for(i in 1:length(t))
        {
            s <- t[i]
            R <- lower.part(L) + s*lower.part(X)
            diag(R) <- diag(L) * exp(s*diag(X)/diag(L))
            res[,,i] <- R
        }
        
        return(res)
    }
    else
    {
        R <- lower.part(L) + t*lower.part(X)
        diag(R) <- diag(L) * exp(t*diag(X)/diag(L))
        return(R)
    }
    
}

#' log-cholesky distance on the Cholesky factor space
#' @keywords internal
geo.dist.chol <- function(L1,L2)
{
    if(is.matrix(L1))
    {
        R <- frobenius(lower.part(L1)-lower.part(L2))
        R <- R + sum((log(diag(L1))-log(diag(L2)))^2)
        return(sqrt(R))
    }
    else # vector form
    {
        n <- length(L1)
        #d <- get.matrix.dim(n)
        d <- round(-0.5+sqrt(2*length(L1)+0.25))
        R <- sum(((L1[(d+1):n])-(L2[(d+1):n]))^2)
        R <- R + sum((log((L1[1:d]))-log((L2[1:d])))^2)
        return(sqrt(R))
    }
    
}

#' log-cholesky log map on the Cholesky factor space
#' @keywords internal
rie.log.chol <- function(L1,L2)
{
    R <- lower.part(L2) - lower.part(L1)
    diag(R) <- diag(L1) * (log(diag(L2))-log(diag(L1)))
    return(R)
}

#' log-cholesky exponential map on the Cholesky factor space
#' @keywords internal
rie.exp.chol <- function(L,X)
{
    return(geodesic.chol(L,X,1))
}

#' The Cholesky factor (lower triangular) of an SPD
#' @param S a SPD
#' @export
ch <- function(S)
{
    return(t(chol(S)))
}

#' Map a Cholesky factor into its corresponding SPD
#' @param L a lower triangular matrix
#' @return \code{L} multiplied by the transpose of \code{L}
#' @keywords internal
ch.to.spd <- function(L)
{
    R <- L %*% t(L)
    return(make.sym(R))
}

#' the differential of ch at S
#' @param S spd
#' @param W symmetric
#' @param L L=ch(S) if it is provided for speedup
#' @keywords internal
d.ch <- function(S,W,L=NULL)
{
    if(is.null(L)) L <- ch(S)
    iL <- solve(L)
    return( L %*%  half(iL %*% W %*% t(iL)))
}

#' the differential of ch.to.spd at L
#' @keywords internal
d.ch.to.spd <- function(L,X)
{
    R <- L %*% t(X) + X %*% t(L)
    return(make.sym(R))
}

#' Log-Cholesky metric of SPD at p
#' @param p spd, the base point
#' @param u symmetric, a tangent vector at p
#' @param v symmetric, a tangent vector at p
#' @param ... L=ch(p)
#' @keywords internal
#' @export
rie.metric.spd_LogCholesky <- function(mfd,p,u,v,...)
{
    opt.param <- list(...)
    if('L' %in% names(opt.param)) L <- opt.param$L
    else L <- ch(p)
    return(metric.chol(L,d.ch(p,u,L),d.ch(p,v,L)))
}

#' Log-Cholesky geodesic of SPD at S
#' @param p spd
#' @param u symmetric
#' @param t time
#' @param ... L=ch(S)
#' @keywords internal
#' @export
geodesic.spd_LogCholesky <- function(mfd,p,u,t,...)
{
    opt.param <- list(...)
    if('L' %in% names(opt.param)) L <- opt.param$L
    else L <- ch(p)
    
    gL <- geodesic.chol(L,d.ch(p,u,L),t)
    if(length(t) > 1)
    {
        R <- array(0,dim(gL))
        for(i in 1:length(t))
            R[,,i] <- make.sym(gL[,,i] %*% base::t(gL[,,i]))
    }
    else
    {
        R <- make.sym(gL %*% base::t(gL))
    }
    
    return(R)

}

#' Compute the geodesic distance of two SPD in Log-Cholesky metric
#' @param p a SPD
#' @param q a SPD
#' @param ... L1=ch(S1) and L2=ch(S2)
#' @keywords internal
#' @export
geo.dist.spd_LogCholesky <- function(mfd,p,q,...)
{
        opt.param <- list(...)
        if(is.null(opt.param) || !('L1' %in% names(opt.param))) L1 <- ch(p)
        else L1 <- opt.param$L1
        if(is.null(opt.param) || !('L2' %in% names(opt.param))) L2 <- ch(q)
        else L2 <- opt.param$L2
        return(geo.dist.chol(L1,L2))
}

#' Riemannian exponential map
#' @param p the foot point of the exp map
#' @param v a tangent vector at \code{S}
#' @param ... L=ch(S)
#' @keywords internal
#' @export
rie.exp.spd_LogCholesky <- function(mfd,p,v,...)
{
    opt.param <- list(...)
    if('L' %in% names(opt.param)) L <- opt.param$L
    else L <- ch(p)
    
    eL <- rie.exp.chol(L,d.ch(p,v,L))
    return(make.sym(eL %*% t(eL)))
}

#' Riemannian logarithmic map
#' @param p the foot point of the log map
#' @param q the point to be mapped
#' @param ... L1=ch(p) and L2=ch(q)
#' @keywords internal
#' @export
rie.log.spd_LogCholesky <- function(mfd,p,q,...)
{
        opt.param <- list(...)
   
        if(is.null(opt.param) || !('L1' %in% names(opt.param))) L1 <- ch(p)
        else L1 <- opt.param$L1
        if(is.null(opt.param) || !('L2' %in% names(opt.param))) L2 <- ch(q)
        else L2 <- opt.param$L2
        
        X <- rie.log.chol(L1,L2)
        W <- d.ch.to.spd(L1,X)
        return(make.sym(W))
}

#' Parallel transport of a tangent vector from a point to another
#' @param p a SPD
#' @param q a SPD
#' @param v a tangent vector at \code{S1}
#' @keywords internal
#' @export
parallel.transport.spd_LogCholesky <- function(mfd,p,q,v,...)
{
    L <- ch(p)
    K <- ch(q)
    X <- d.ch(p,v,L)
    Y <- lower.part(X,T) + diag(as.vector(diag(K)*diag(L)^(-1)*diag(X)),nrow=mfd$dim[1],ncol=mfd$dim[2])
    return(d.ch.to.spd(K,Y))
}




#' Generate a set of normally distributed random matrices that represent tangent vectors at some point. 
#' @param mfd an object created by \code{matrix.manifold}
#' @param n sample size
#' @param sig the standard deviation of the normal distribution
#' @param drop whether return the result as a matrix when \code{n=1}
#' @return an \code{M*N*n} array of \code{n} matrices, where \code{M*N} is the dimensions of matrices
#' @keywords internal
#' @export
rtvecor.spd_LogCholesky <- function(mfd,n=1,sig=1,drop=T)
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
#' @keywords internal
#' @export
rmatrix.spd_LogCholesky <- function(mfd,n=1,mu=NULL,sig=1,drop=T)
{
    
    if(is.null(mu)) mu <- diag(rep(1,mfd$dim[1]))
    stopifnot(is.spd(mu))
    
    S <- rtvecor(mfd,n,sig,drop=F)
    
    R <- array(0,c(mfd$dim[1],mfd$dim[1],n))
    for(i in 1:n)
    {
        R[,,i] <- rie.exp(mfd,mu,S[,,i])
    }
    
    if(n==1 && drop) return(R[,,i])
    else return(R)
}

#' Frechet mean of matrices
#' @param mfd an object created by \code{matrix.manifold}
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