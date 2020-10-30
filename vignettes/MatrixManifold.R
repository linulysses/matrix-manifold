## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
library(MatrixManifold)
mfd <- matrix.manifold('spd','LogCholesky',c(3,3))
S <- rmatrix(mfd,n=10,mu=rmatrix(mfd,n=1,diag(3)))
mu <- frechet.mean(mfd,S)
mu

# P <- rmatrix(mfd)
# Q <- rmatrix(mfd)
# geo.dist(mfd,P,Q)
# 
# t <- runif(1)
# V <- rtvecor(mfd)
# geodesic(mfd,P,V,t)
# 
# V <- rie.log(mfd,P,Q)
# Y <- rie.exp(mfd,P,V)
# Y

