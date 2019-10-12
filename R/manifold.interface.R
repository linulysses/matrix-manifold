#' Create an object that represent a Riemannian structure on a matrix space
#' @param dim the dimensions of the matrices, must be an positive integer or a vector of two positive integers
#' @param manifold matrix manifold, one of 'lt', 'spd','stiefel' and 'sym'
#' @param metric Riemannian metric, 'LogCholesky','AffineInvariant','LogEuclidean' for spd
#' @return an object containing essential information about the defined matrix manifold
#' @export
create.matrix.manifold <- function(manifold,metric,dim)
{
    stopifnot(manifold %in% c('spd')) #,'stiefel','sym'))
    stopifnot(metric %in% c('LogCholesky'))#,'AffineInvariant','LogEuclidean'))
    class.name <- paste0(manifold,'_',metric)
    if(length(dim) == 1) dim <- rep(dim,2)
    mfd <- structure(list(dim=dim,manifold=manifold,metric=metric),class=class.name)
    return(mfd)
}

#' Compute the Riemannian exponential map
#' @param mfd the manifold object created by \code{create.matrix.manifold}
#' @param p the foot point of the exponential map
#' @param v a tangent vector at \code{p}
#' @param ... other parameters (primarily for computation speedup)
#' @return a matrix that is the exponential map of \code{v}
#' @export
rie.exp <- function(mfd,p,v,...)
{
    UseMethod("rie.exp")
}

#' Compute the Riemannian logarithm map
#' @param mfd the manifold object created by \code{create.matrix.manifold}
#' @param p the foot point of the logarithm map
#' @param q a point on the manifold
#' @param ... other parameters (primarily for computation speedup)
#' @return a matrix that is the logarithm map of \code{q}
#' @export
rie.log <- function(mfd,p,q,...)
{
    UseMethod('rie.log')
}

#' Compute the Riemannian metric 
#' @param mfd the manifold object created by \code{create.matrix.manifold}
#' @param p the foot point of the metric
#' @param u a tangent vector at \code{p}
#' @param v a tangent vector at \code{q}
#' @param ... other parameters (primarily for computation speedup)
#' @return a number that is the Riemannian metric of \code{u} and \code{v}
#' @export
rie.metric <- function(mfd,p,u,v,...)
{
    UseMethod("rie.metric")
}

#' Compute the geodesic starting at some point on a manifold
#' @param mfd the manifold object created by \code{create.matrix.manifold}
#' @param p the starting point of the geodesic
#' @param u a tangent vector at \code{p}, representing the direction of the geodesic
#' @param t a vector of nonnegative real numbers
#' @param ... other parameters (primarily for computation speedup)
#' @return an array of matrices which are the geodesic evaluated at \code{t}
#' @export
geodesic <- function(mfd,p,u,t,...)
{
    UseMethod('geodesic')
}

#' Compute the geodesic distance
#' @param mfd the manifold object created by \code{create.matrix.manifold}
#' @param p a point on the manifold
#' @param q a point on the manifold
#' @param ... other parameters (primarily for computation speedup)
#' @return the geodesic distance between \code{p} and \code{q}
#' @export
geo.dist <- function(mfd,p,q,...)
{
    UseMethod('geo.dist')
}

#' Parallel transport of a tangent vector from one point to another along the geodesic between the two points
#' @param mfd the manifold object created by \code{create.matrix.manifold}
#' @param p a point on the manifold
#' @param q a point on the manifold
#' @param v a tangent vector at \code{p}
#' @param ... other parameters (primarily for computation speedup)
#' @return a tangent vector at \code{q}
#' @export
parallel.transport <- function(mfd,p,q,v,...) 
{
    UseMethod('parallel.transport')
}

#' Generate a set of random matrices that represent tangent vectors at some point. 
#' @param mfd an object created by \code{create.matrix.manifold}
#' @param n sample size
#' @param sig the standard deviation of the normal distribution
#' @return an \code{M*N*n} array of \code{n} matrices, where \code{M*N} is the dimensions of matrices
#' @details The generated samples have expectation zero and follow a isotropic D-dimensional normal distribution with isotropic variance \code{sig}, where D is the intrinsic dimension of the matrix manifold
#' @export
gen.tangent.vectors <- function(mfd,n=1,sig=1,drop=T)
{
    UseMethod('gen.tangent.vectors')
}

#' Generate a set of random matrices on a matrix manifold
#' @param mfd an object created by \code{create.matrix.manifold}
#' @param n sample size
#' @param mu the Frechet mean. If \code{NULL} is given, then it is the identity element, e.g., for SPD, it is the identity matrix
#' @param sig the standard deviation of the normal distribution
#' @return an \code{M*N*n} array of \code{n} matrices, where \code{M*N} is the dimensions of matrices
#' @details The generated samples have Frechet mean \code{mu}. The logarithmic maps of these samples at \code{mu} follow a isotropic D-dimensional normal distribution with isotropic variance \code{sig}, where D is the intrinsic dimension of the matrix manifold
#' @export
gen.matrices <- function(mfd,n=1,mu=NULL,sig=1,drop=T)
{
    UseMethod('gen.matrices')
}

#' Frechet mean of matrices
#' @param mfd an object created by \code{create.matrix.manifold}
#' @param S an \code{M*N*n} array of matrices or a list of \code{n} matrices, where \code{n} is the number of matrices
#' @return the Frechet mean of the matrices in \code{S}
#' @export
frechet.mean <- function(mfd,S)
{
    UseMethod('frechet.mean')
}