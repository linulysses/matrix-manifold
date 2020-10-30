context("check routines in triangular.matrix.R")

test_that("check routines in triangular.matrix.R", {

    eps <- 1e-2
    
    for(d in c(1,3))
    {
        mfd <- matrix.manifold('ltpd','LogCholesky',dim=c(d,d))
        
        A <- rmat(d,d)
        P <- make.sym(A %*% t(A))
        L <- ch(P)
        expect_true(is.lower(L))
        expect_true(is.ltpd(L))
        
        W <- rtvecor.ltpd_LogCholesky(mfd,1)
        V <- rtvecor.ltpd_LogCholesky(mfd,1)
        m <- rie.metric.ltpd_LogCholesky(mfd,L,W,V)
        
        expect_true(abs(m-metric.chol(L,W,V)) < eps)
        
        W <- runif(1) * W / sqrt(rie.metric.ltpd_LogCholesky(mfd,P,W,W))
        S <- rie.exp.ltpd_LogCholesky(mfd,P,W)
        V <- rie.log.ltpd_LogCholesky(mfd,P,S)
        expect_true(frobenius.norm(W-V) < eps)
        
        t <- runif(1)
        W <- W / sqrt(rie.metric.ltpd_LogCholesky(mfd,P,W,W))
        S <- geodesic.ltpd_LogCholesky(mfd,P,W,t)
        expect_true(abs(geo.dist.ltpd_LogCholesky(mfd,P,S)-t) < eps)
        
        
        W <- W / sqrt(rie.metric.ltpd_LogCholesky(mfd,P,W,W))
        t <- runif(2)
        Qs <- geodesic.ltpd_LogCholesky(mfd,P,W,t)
        Q <- as.matrix(Qs[,,1])
        X <- parallel.transport.ltpd_LogCholesky(mfd,P,Q,V)
        Y <- parallel.transport.ltpd_LogCholesky(mfd,P,Q,W)
        expect_true(abs(rie.metric.ltpd_LogCholesky(mfd,P,W,V)-
                            rie.metric.ltpd_LogCholesky(mfd,Q,X,Y)) < eps)
        
        mfd <- matrix.manifold('ltpd','LogCholesky',c(d,d))
        
        for(n in c(1,10))
        {
            mu <- diag(rep(1,d)) 
            
            S <- rmatrix.ltpd_LogCholesky(mfd,n,drop=F,mu=mu)
            expect_true(all(sapply(1:n, function(i){
                is.ltpd(as.matrix(S[,,i]))
            })))
            V <- rtvecor.ltpd_LogCholesky(mfd,n=n,sig=0.1,drop=F)
            expect_true(all(sapply(1:n, function(i){
                is.lower(as.matrix(V[,,i]))
            })))
            
            mu <-rmatrix.ltpd_LogCholesky(mfd,n=1)
            V <- center.matrices(V)
            Q <- array(0,c(d,d,n))
            for(i in 1:n)
            {
                Q[,,i] <- rie.exp.ltpd_LogCholesky(mfd,mu,as.matrix(V[,,i]))
            }
            Q.mu <- frechet.mean.ltpd_LogCholesky(mfd,Q)
            
            expect_true(frobenius.norm(mu-Q.mu) < eps)
        }
        
    }
    
    
    
})