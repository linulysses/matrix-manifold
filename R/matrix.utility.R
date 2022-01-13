#' matrix inverse
#' @keywords internal
inv <- function(S)
{
    solve(S)
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


#' Generate a random matrix
#' @keywords internal
rmat <- function(m,n)
{
    matrix(runif(m*n),m,n)
}

#' Check whether a matrix is a symmetric positive definite matrix 
#' @param A a matrix
#' @return TRUE or FALSE depending on whether \code{A} is an SPD
#' @keywords internal
is.spd <- function(A)
{
    if(is.vec(A) && length(A)==1) return(A>0)
    return(is.sym(A) && matrixcalc::is.positive.definite(A))
}

#' Check whether a matrix is a symmetric matrix 
#' @param A a matrix
#' @return TRUE or FALSE depending on whether \code{A} is symmetric
#' @keywords internal
is.sym <- function(A)
{
    A <- as.matrix(A)
    return(matrixcalc::is.symmetric.matrix(A))
}

#' Check whether an object is a vector
#' @keywords internal
is.vec <- function(v)
{
    is.vector(v) || (is.matrix(v) && any(dim(v)==1))
}

#' Check whether a matrix ia a lower triangular matrix
#' @keywords internal
is.lower <- function(L)
{
    all(strict.upper.part(L)==0)
}


#' Check whether a matrix ia a LTPD
#' @keywords internal
is.ltpd <- function(L)
{
    all(strict.upper.part(L)==0) & all(diag(L)>0)
}

#' Make a matrix into a symmetric matrix to correct potential asymmetry due to finite-precision operations
#' @keywords internal 
make.sym <- function(S)
{
    return((S+t(S))/2)
}


#' Frobenius inner product of two matrices
#' @param A a matrix
#' @param B a matrix with the same dimensions of \code{A}. If \code{B} is not provided, it is set to \code{A}
#' @return the Frobenius inner product of \code{A} and \code{B}
#' @keywords internal
frobenius <- function(A,B=NULL)
{
    if(is.null(B)) B <- A
    return(matrixcalc::frobenius.prod(A,B))
}

#' Frobenius norm of a matrix
#' @param A a matrix
#' @return the Frobenius norm of \code{A}
#' @keywords internal
frobenius.norm <- function(A)
{
    return(matrixcalc::frobenius.norm(A))
}

#' return the lower triangular part of A
#' @param strict whether a strict triangular matrix is returned, i.e. whether the diagonal elements are EXCLUDED (strict) or not (non-strict)
#' @keywords internal
lower.part <- function(A,strict=T)
{
        A <- matrixcalc::lower.triangle(A)
        if(strict) diag(A) <- 0
        return(A)

    
}

#' The strict lower triangular part of a matrix
#' @keywords internal
strict.lower.part <- function(A)
{
    return(lower.part(A,T))
}

#' The diagonal matrix of the diagonal part of a matrix
#' @keywords internal
diag.part <- function(A)
{
    return(diag(diag(A)))
}

#' The upper triangular part of a matrix
#' @param strict whether a strict triangular matrix is returned, i.e. whether the diagonal elements are EXCLUDED (strict) or not (non-strict)
#' @keywords internal
upper.part <- function(A,strict=T)
{
    A <- matrixcalc::upper.triangle(A)
    if(strict) diag(A) <- 0
    return(A)
}

#' The strict upper triangular part of a matrix
#' @keywords internal
strict.upper.part <- function(A)
{
    return(upper.part(A,T))
}

#' the half operation of A
#' return a matrix B of the same dimension of A, such that B(i,j)=A(i,j) if i>j, B(i,i)=A(i,i)/2, and B(i,j)=0 if i<j
#' @keywords internal
half <- function(A)
{
    R <- lower.part(A)
    diag(R) <- diag(A) / 2
    return(R)
}

#' Matrix exponential of a diagonal matrix
#' @keywords internal
expm.diag <- function(D)
{
    D <- as.matrix(D)
    diag(D) <- exp(diag(D))
    return(D)
}

#' Matrix log of a diagonal matrix
#' @keywords internal
logm.diag <- function(D)
{
    diag(D) <- log(diag(D))
    return(D)
}

#' Inverse of a diagonal matrix
#' @keywords internal
inv.diag <- function(D)
{
    diag(D) <- 1/diag(D)
    return(D)
}

#' Multiplication of two diagonal matrices
#' @keywords internal
mul.diag <- function(D1,D2)
{
    diag(D1) <- diag(D1) * diag(D2)
    return(D1)
}

#' Trace of a matrix
#' @param S square matrix
#' @keywords internal
trace <- function(S)
{
    sum(diag(S))
}

#' 
#' 
#' #' Euler 3-by-3 rotation matrix with the given angles 
#' #' @keywords internal
#' euler.rot <- function(phi,theta,psi)
#' {
#'     Rot <- matrix(0,3,3)
#'     cosphi = cos(phi)
#'     cospsi = cos(psi)
#'     costheta = cos(theta)
#'     sinphi = sin(phi)
#'     sinpsi = sin(psi)
#'     sintheta = sin(theta)
#'     
#'     Rot[1,1] = cosphi*costheta*cospsi - sinphi*sinpsi
#'     Rot[2,1] = sinphi*costheta*cospsi + cosphi*sinpsi
#'     Rot[3,1] = -sintheta*cospsi
#'     
#'     Rot[1,2] = -cosphi*costheta*sinpsi - sinphi*cospsi
#'     Rot[2,2] = -sinphi*costheta*sinpsi + cosphi*cospsi
#'     Rot[3,2] = sintheta*sinpsi
#'     
#'     Rot[1,3] = cosphi*sintheta
#'     Rot[2,3] = sinphi*sintheta
#'     Rot[3,3] = costheta
#'     return(Rot)
#' }
#' 
#' 


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

