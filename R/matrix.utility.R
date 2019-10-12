#' Generate a random matrix
#' @keywords internal
rmatrix <- function(m,n)
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
