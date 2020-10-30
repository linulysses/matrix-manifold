context("check routines in symmetric.matrix.R")

test_that("check routines in symmetric.matrix.R", {
    
    eps <- 1e-12
    
    for(d in c(1,3))
    {
        for(n in c(1,10))
        {
            L <- make.sym(rmat(d,d))
            v <- sym.to.vec(L,drop=T)
            U <- vec.to.sym(v,drop=T)
            expect_true(all(L==U))
            
            S <- rsym(d,n,drop=F)
            V <- sym.to.vec(S,drop=F)
            expect_true(all(dim(V)==c(round(d*(d+1)/2),n)))
            R <- vec.to.sym(V,drop=F)
            expect_true(all(S==R))
            
            expect_true(is.sym.in.vec.form(V))
            expect_false(is.sym.in.vec.form(S))
            
            expect_true(is.sym.in.matrix.form(S))
            expect_false(is.sym.in.matrix.form(V))
            
            S <- make.sym(rmat(d,d))
            expect_true(is.sym.in.matrix.form(S))
            
            S1 <- sym.to.matrix.form(S)
            expect_true(all(S==S1))
            
            S <- rsym(d,n,drop=F)
            S1 <- sym.to.matrix.form(S)
            expect_true(all(S==S1))
            
            V <- sym.to.vec.form(S)
            S1 <- sym.to.matrix.form(V)
            expect_true(all(S==S1))
            
            m <- count.sym.matrix(S)
            expect_true(m==n)
        }
        
        mfd <- matrix.manifold('sym','Frobenius',dim=c(d,d))
        
        A <- rmat(d,d)
        P <- make.sym(A %*% t(A))
        W <- rsym(d,1)
        V <- rsym(d,1)
        m <- rie.metric.sym_Frobenius(mfd,P,W,V)
        expect_true(abs(m-sum(W*V)) < eps)
        
        W <- runif(1) * W / sqrt(rie.metric.sym_Frobenius(mfd,P,W,W))
        S <- rie.exp.sym_Frobenius(mfd,P,W)
        V <- rie.log.sym_Frobenius(mfd,P,S)
        expect_true(frobenius.norm(W-V) < eps)
        
        t <- runif(1)
        W <- W / sqrt(rie.metric.sym_Frobenius(mfd,P,W,W))
        S <- geodesic.sym_Frobenius(mfd,P,W,t)
        expect_true(abs(geo.dist.sym_Frobenius(mfd,P,S)-t) < eps)
        
        
        A <- rmat(d,d)
        P <- make.sym(A %*% t(A))
        W <- rsym(d,1)
        W <- W / sqrt(rie.metric.sym_Frobenius(mfd,P,W,W))
        V <- rsym(d,1)
        t <- runif(2)
        Qs <- geodesic.sym_Frobenius(mfd,P,W,t)
        Q <- Qs[,,1]
        X <- parallel.transport.sym_Frobenius(mfd,P,Q,V)
        Y <- parallel.transport.sym_Frobenius(mfd,P,Q,W)
        expect_true(abs(rie.metric.sym_Frobenius(mfd,P,W,V)-
                            rie.metric.sym_Frobenius(mfd,Q,X,Y)) < eps)
        
    }
    
    
})