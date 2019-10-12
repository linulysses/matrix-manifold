context("check routines in log.cholesky.R")

test_that("check routines in log.cholesky.R", {
    
    eps <- 1e-2
    
    for(d in c(1,3))
    {
        A <- rmatrix(d,d)
        P <- make.sym(A %*% t(A))
        L <- ch(P)
        expect_true(is.lower(L))
        expect_true(frobenius.norm(P-L%*%t(L)) < eps)
        expect_true(frobenius.norm(P-ch.to.spd(L)) < eps)
        
        X <- lower.part(A,strict=F)
        W <- d.ch.to.spd(L,X)
        Y <- d.ch(ch.to.spd(L),W)
        expect_true(frobenius.norm(X-Y) < eps)
        
        
        R <- geodesic.chol(L,X,1)
        expect_true(frobenius.norm(rie.exp.chol(L,X)-R) < eps)
        expect_true(frobenius.norm(X-rie.log.chol(L,R)) < eps)
        
        t <- runif(1)
        X <- X / sqrt(metric.chol(L,X,X))
        R <- geodesic.chol(L,X,t)
        dist <- geo.dist.chol(L,R)
        expect_true(abs(dist-t) < eps)
        
        
        mfd <- create.matrix.manifold('spd','LogCholesky',dim=c(d,d))
        
        A <- rmatrix(d,d)
        P <- make.sym(A %*% t(A))
        W <- gen.sym(d,1)
        V <- gen.sym(d,1)
        m <- rie.metric.spd_LogCholesky(mfd,P,W,V)
        L <- ch(P)
        X <- d.ch(P,W)
        Y <- d.ch(P,V)
        expect_true(abs(m-metric.chol(L,X,Y)) < eps)
        
        W <- runif(1) * W / sqrt(rie.metric.spd_LogCholesky(mfd,P,W,W))
        S <- rie.exp.spd_LogCholesky(mfd,P,W)
        V <- rie.log.spd_LogCholesky(mfd,P,S)
        expect_true(frobenius.norm(W-V) < eps)
        
        t <- runif(1)
        W <- W / sqrt(rie.metric.spd_LogCholesky(mfd,P,W,W))
        S <- geodesic.spd_LogCholesky(mfd,P,W,t)
        expect_true(abs(geo.dist.spd_LogCholesky(mfd,P,S)-t) < eps)
        
        
        A <- rmatrix(d,d)
        P <- make.sym(A %*% t(A))
        W <- gen.sym(d,1)
        W <- W / sqrt(rie.metric.spd_LogCholesky(mfd,P,W,W))
        V <- gen.sym(d,1)
        t <- runif(1)
        Q <- geodesic.spd_LogCholesky(mfd,P,W,t)
        X <- parallel.transport.spd_LogCholesky(mfd,P,Q,V)
        Y <- parallel.transport.spd_LogCholesky(mfd,P,Q,W)
        expect_true(abs(rie.metric.spd_LogCholesky(mfd,P,W,V)-
                            rie.metric.spd_LogCholesky(mfd,Q,X,Y)) < eps)
        
    }
    
    
})