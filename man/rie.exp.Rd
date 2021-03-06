% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manifold.interface.R
\name{rie.exp}
\alias{rie.exp}
\title{Compute the Riemannian exponential map}
\usage{
rie.exp(mfd, p, v, ...)
}
\arguments{
\item{mfd}{the manifold object created by \code{matrix.manifold}}

\item{p}{the foot point of the exponential map}

\item{v}{a tangent vector at \code{p}}

\item{...}{other parameters (primarily for computation speedup)}
}
\value{
a matrix that is the exponential map of \code{v}
}
\description{
Compute the Riemannian exponential map
}
