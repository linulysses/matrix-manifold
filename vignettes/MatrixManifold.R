## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
library(MatrixManifold)
mfd <- create.matrix.manifold('spd','LogCholesky',c(3,3))
S <- gen.matrices(mfd,n=10,mu=gen.matrices(mfd,n=1,diag(3)))
mu <- frechet.mean(mfd,S)
mu

P <- gen.matrices(mfd)
Q <- gen.matrices(mfd)
geo.dist(mfd,P,Q)

t <- runif(1)
V <- gen.tangent.vectors(mfd)
geodesic(mfd,P,V,t)

V <- rie.log(mfd,P,Q)
Y <- rie.exp(mfd,P,V)
Y

