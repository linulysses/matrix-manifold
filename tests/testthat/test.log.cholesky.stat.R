context("check routines in log.cholesky.stat.R")

test_that("check routines in log.cholesky.stat.R", {
    
    eps <- 1e-12
    
    for(d in c(1,3))
    {
        mfd <- create.matrix.manifold('spd','LogCholesky',c(d,d))
        
        for(n in c(1,10))
        {
            mu <- diag(rep(1,d)) #gen.matrices.spd_LogCholesky(mfd,n=1,sig=0.1,drop=T)
            
            S <- gen.sym(d,n,drop=F)
            expect_true(all(dim(S)==c(mfd$dim,n)))
            
            S <- gen.matrices.spd_LogCholesky(mfd,n,drop=F,mu=mu)
            expect_true(all(sapply(1:n, function(i){
                is.spd(S[,,i])
            })))
            V <- gen.tangent.vectors.spd_LogCholesky(mfd,n=n,sig=0.1,drop=F)
            expect_true(all(sapply(1:n, function(i){
                is.sym(V[,,i])
            })))
            
            mu <- gen.matrices.spd_LogCholesky(mfd,n=1)
            V <- center.matrices(V)
            Q <- array(0,c(d,d,n))
            for(i in 1:n)
            {
                Q[,,i] <- rie.exp.spd_LogCholesky(mfd,mu,V[,,i])
            }
            Q.mu <- frechet.mean.spd_LogCholesky(mfd,Q)
            
            expect_true(frobenius.norm(mu-Q.mu) < eps)
        }
        
        
    }
    
    
})