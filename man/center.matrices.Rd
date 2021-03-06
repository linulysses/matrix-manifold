% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix.utility.R
\name{center.matrices}
\alias{center.matrices}
\title{#' Euler 3-by-3 rotation matrix with the given angles 
#' @keywords internal
euler.rot <- function(phi,theta,psi)
{
    Rot <- matrix(0,3,3)
    cosphi = cos(phi)
    cospsi = cos(psi)
    costheta = cos(theta)
    sinphi = sin(phi)
    sinpsi = sin(psi)
    sintheta = sin(theta)
    
    Rot[1,1] = cosphi*costheta*cospsi - sinphi*sinpsi
    Rot[2,1] = sinphi*costheta*cospsi + cosphi*sinpsi
    Rot[3,1] = -sintheta*cospsi
    
    Rot[1,2] = -cosphi*costheta*sinpsi - sinphi*cospsi
    Rot[2,2] = -sinphi*costheta*sinpsi + cosphi*cospsi
    Rot[3,2] = sintheta*sinpsi
    
    Rot[1,3] = cosphi*sintheta
    Rot[2,3] = sinphi*sintheta
    Rot[3,3] = costheta
    return(Rot)
}}
\usage{
center.matrices(S, mu = NULL)
}
\arguments{
\item{S}{an array of matrices where the last dimension corresponds to sample size}

\item{mu}{the center}
}
\description{
Center a sample of matrices in Euclidean space
}
\keyword{internal}
