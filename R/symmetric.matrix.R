
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

#' Transform a lower triangular matrix into a vector representation, diagonal by diagonal, starting from the main diagonal
#' @keywords internal
sym.to.vec.atomic <- function(S)
{
    if(is.vec(S) && length(S)==1) return(S)
    
    d <- dim(S)[1]
    v <- rep(0,d*(d+1)/2)
    v[1:d] <- diag(S)
    p <- d
    for(i in 2:d)
        for(j in 1:(d-i+1))
        {
            p <- p + 1
            v[p] <- S[i+j-1,j]
        }
    attr(v,'format') <- 'sym.in.vec'
    return(v)         
}

#' Transform a lower triangular matrix into a vector representation, diagonal by diagonal, starting from the main diagonal
#' @keywords internal
sym.to.vec <- function(S,drop=F)
{
    if(is.matrix(S)) #return(sym.to.vec.atomic(S))
        S <- array(S,dim=c(dim(S),1))
    
    if(is.array(S))
    {
        n <- count.sym.matrix(S)
        d <- dim(S)[1]
        V <- matrix(0,d*(d+1)/2,n)
        for(i in 1:n)
        {
            V[,i] <- sym.to.vec.atomic(S[,,i])
        }
        if(drop && n==1) V <- V[,1] 
        attr(V,'format') <- 'sym.in.vec'
        return(V)
    }
    else if(is.list(S))
    {
        lapply(S,sym.to.vec.atomic)
    }
}


#' Transform a vector created by \code{sym.to.vec.atomic} or \code{sym.to.vec} into a symmetric matric
#' @keywords internal
vec.to.sym.atomic <- function(v)
{
    d <- round(-0.5+sqrt(2*length(v)+0.25))
    S <- matrix(0,d,d)
    p <- 0
    for(i in 1:d)
        for(j in 1:(d-i+1))
        {
            p <- p + 1
            S[i+j-1,j] <- v[p]
            S[j,i+j-1] <- v[p]
        }
    attr(S,'format') <- 'sym.in.matrix'
    return(S)         
}

#' recover a (set of) lower triangular matrix from its vector representation
#' @keywords internal
vec.to.sym <- function(v,drop=F)
{
    if(is.vec(v))
    {
        v <- as.matrix(v)
        attr(v,'format') <- 'sym.in.vec'
    }
        
    
    if(is.matrix(v))
    {
        stopifnot(attr(v,'format') == 'sym.in.vec')
        
        n <- ncol(v)
        d <- round(-0.5+sqrt(2*length(v[,1])+0.25))
        S <- array(0,c(d,d,n))
        for(i in 1:n)
        {
            S[,,i] <- vec.to.sym.atomic(v[,i])
        }
        if(drop && n == 1) return(S[,,1])
        else return(S)
    }
    else if(is.list(v))
    {
        lapply(v,vec.to.sym.atomic)
    }
}

#' Check whether a symmetric matrix is represented by its vector form
#' @keywords internal
is.sym.in.vec.form <- function(S)
{
    if(is.list(S))
    {
        return(is.sym.in.vec.form(S[[1]]))
    }
    else if(is.array(S) || is.matrix(S))
    {
        return(attr(S,'format') == 'sym.in.vec')
    }
    else if(is.vec(S) && length(S) == 1) return(TRUE)
    else return(FALSE)
}

#' Check whether a matrix is represented by its vector form
#' @keywords internal
is.sym.in.matrix.form <- function(S)
{
    if(is.list(S))
    {
        return(is.sym.in.matrix.form(S[[1]]))
    }
    else if(!is.null(attr(S,'format')))
    {
        return(attr(S,'format') == 'sym.in.matrix')
    }
    else if(is.matrix(S))
    {
        return(is.sym(S))
    }
    else if(is.array(S))
    {
        return(is.sym.in.matrix.form(S[,,1]))
    }
    else return(FALSE)
}

#' transform a matrix into its matrix form
#' @keywords internal
sym.to.matrix.form <- function(S)
{
    if(is.sym.in.matrix.form(S)) return(S)
    else if(is.sym.in.vec.form(S))
    {
        return(vec.to.sym(S))
    }
    else stop('unrecoginzed format of S')
}

#' transform a matrix into its vector form
#' @keywords internal
sym.to.vec.form <- function(S)
{
    if(is.sym.in.vec.form(S)) return(S)
    else if(is.sym.in.matrix.form(S))
    {
        return(sym.to.vec(S))
    }
    else stop('unrecoginzed format of S')
}

#' get the matrix dimensions 
#' @keywords internal
get.sym.matrix.dim <- function(v)
{
    if(is.list(v)) return(get.sym.matrix.dim(v[[1]]))
    
    if(is.sym.in.vec.form(v))
    {
        if(is.matrix(v)) n <- nrow(v)
        else if(is.vec(v)) n <- length(v)
        else stop('unrecognized v')
        d <- round(-0.5+sqrt(2*n+0.25))
        return(d)
    }
    else if(is.sym.in.matrix.form(v))
    {
        if(is.array(v)) v <- v[,,1]
        return(nrow(v))
    }
    if(is.matrix(v)) return(nrow(v))
    else if(is.vec(v))
    {
        n <- length(v)
        d <- round(-0.5+sqrt(2*n+0.25))
        return(d)
    }
    else stop('unrecognized v')

}

#' count the number of matrices represented by an array or a list
#' @keywords internal
count.sym.matrix <- function(S)
{
    format <- attr(S,'format')
    if(!is.null(format) && format == 'sym.in.vec') return(ncol(S))
    else if(!is.null(format) && format == 'sym.in.matrix')
    {
        if(is.matrix(S)) return(1)
        else if(is.array(S)) return(dim(S)[3])
    }
    else if(is.array(S)) return(dim(S)[3])
    else if(is.matrix(S) && is.sym(S)) return(1)
    else if(is.list(S)) return(length(S))
    else stop('unrecognized format of S')
}



#' Frobenius metric of symmetric matrices at S
#' @param S symmetric, the base point
#' @param W symmetric, a tangent vector at S
#' @param V symmetric, a tangent vector at S
#' @param opt.param optional parameters, primarily for speedup. Not used here
#' @keywords internal
#' @export
rie.metric.sym_Frobenius <- function(mfd,S,W,V,opt.param=NULL)
{
    return(sum(W*V))
}

#' Frobenius exponential map of symmetric matrices at S
#' @param S the foot point of the exponential map
#' @param V a tangent vector at \code{S}
#' @param ... other parameters (primarily for computation speedup)
#' @return a matrix that is the exponential map of \code{V}
#' @export
rie.exp.sym_Frobenius <- function(mfd,S,V,...)
{
    return(S+V)
}

#' Frobenius logarithmic map of symmetric matrices at S
#' @param S the foot point of the log map
#' @param Q a tangent vector at \code{S}
#' @param ... other parameters (primarily for computation speedup)
#' @return a matrix that is the log map of \code{Q}
#' @export
rie.log.sym_Frobenius <- function(mfd,S,Q,...)
{
    return(Q-S)
}

#' Frobenius geodesic of symmetric matrices at S
#' @param S symmetric
#' @param W symmetric
#' @param t time
#' @param opt.param optional parameters, primarily for speedup. Not used here.
#' @keywords internal
#' @export
geodesic.sym_Frobenius <- function(mfd,S,W,t,opt.param=NULL)
{
    if(length(t) == 1) return(S+t*W)
    else
    {
        R <- array(0,c(dim(S),length(t)))
        for(i in 1:length(t))
        {
            R[,,i] <- S + t[i] * W
        }
        return(R)
    }
  
}


#' Compute the geodesic distance of two symmetric in Forbenius metric
#' @param S1 a symmetric
#' @param S2 a symmetric
#' @param opt.param optional parameters, primarily for speedup. 
#' @keywords internal
#' @export
geo.dist.sym_Frobenius <- function(mfd,S1,S2,opt.param=NULL)
{
    return(sqrt(sum((S1-S2)^2)))
}


#' Parallel transport of a tangent vector from a point to another
#' @param S1 symmetric
#' @param S2 symmetric
#' @param W a tangent vector at \code{S1}
#' @keywords internal
#' @export
parallel.transport.sym_Frobenius <- function(mfd,S1,S2,W)
{
    return(W)
}


#' Generate a set of normally distributed random matrices that represent tangent vectors at some point. 
#' @param mfd an object created by \code{create.matrix.manifold}
#' @param n sample size
#' @param sig the standard deviation of the normal distribution
#' @return an \code{M*N*n} array of \code{n} matrices, where \code{M*N} is the dimensions of matrices
#' @keywords internal
#' @export
gen.tangent.vectors.sym_Frobenius <- function(mfd,n=1,sig=1,drop=T)
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
gen.matrices.sym_Frobenius <- function(mfd,n=1,mu=NULL,sig=1,drop=T)
{
    
    if(is.null(mu)) mu <- diag(rep(1,mfd$dim[1]))
    stopifnot(is.sym(mu))
    
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
frechet.mean.sym_Frobenius <- function(mfd,S)
{
    R <- 0
    if(is.list(S))
    {
        for(i in 1:length(S))
        {
            R <- R + S[[i]]
        }
        R <- R / length(S)
        return(make.sym(R))
    }
    else if(is.array(S))
    {
        n <- dim(S)[3]
        d <- dim(S)[1]
        for(i in 1:n)
        {
            R <- R + S[,,i]
        }
        R <- R / n
        return(make.sym(R))
    }
    else if(is.matrix(S)) return(S)
    else stop('S must be an array, a list or a matrix')
    
}