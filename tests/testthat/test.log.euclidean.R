context("check routines in log.euclidean.R")

test_that("check routines in log.euclidean.R", {
    
    eps <- 1e-2
    
    library(expm)
    
    M <- matrix(c(2.23352473,  0.62537226, -0.42983609,
                  0.62537226,  1.34851589, -0.49748003,
                  -0.42983609, -0.49748003,  0.822512), 3,3, byrow=T)
    
    
    dM <- matrix(c(0.77675883, -0.15213052,  1.18733413,
                   -0.15213052,  0.45820972,  0.01693541,
                   1.18733413,  0.01693541, -0.40507469),3,3, byrow=T)
    
    R <- 0
    max.it <- 100
    for(k in 1:max.it)
    {
        S <- 0
        for(j in 0:(k-1)) 
            S <- S + (M %^% (k-j-1)) %*% dM %*% (M %^% j)
        R <- R + S / factorial(k)
    }
    
    expect_true(frobenius.norm(R-d.expm(M,dM)) < eps)
    
    for(d in c(1,3))
    {

        mfd <- matrix.manifold('spd','LogEuclidean',dim=c(d,d))
        
        P <- diag(d)
        W <- rsym(d,1)
        V <- rsym(d,1)
        m <- rie.metric.spd_LogEuclidean(mfd,P,W,V)
        expect_true(abs(m-sum(W*V)) < eps)
        
        A <- rmat(d,d)
        P <- make.sym(A %*% t(A))
        W <- runif(1) * W / sqrt(rie.metric.spd_LogEuclidean(mfd,P,W,W))
        S <- rie.exp.spd_LogEuclidean(mfd,P,W)
        V <- rie.log.spd_LogEuclidean(mfd,P,S)
        expect_true(frobenius.norm(W-V) < eps)
        
        U <- geodesic(mfd,P,V,1)
        expect_true(frobenius.norm(U-S) < eps)
        
        Sf <- geodesic(mfd,P,V,1/2)
        expect_true(abs(geo.dist(mfd,P,Sf)^2-rie.metric(mfd,P,W,W)/4) < eps)
        
        
        A <- rmat(d,d)
        P <- make.sym(A %*% t(A))
        W <- rsym(d,1)
        W <- W / sqrt(rie.metric.spd_LogEuclidean(mfd,P,W,W))
        V <- rsym(d,1)
        A <- rmat(d,d)
        Q <- make.sym(A %*% t(A))
        X <- parallel.transport.spd_LogEuclidean(mfd,P,Q,V)
        Y <- parallel.transport.spd_LogEuclidean(mfd,P,Q,W)
        expect_true(abs(rie.metric.spd_LogEuclidean(mfd,P,W,V)-
                            rie.metric.spd_LogEuclidean(mfd,Q,X,Y)) < eps)
        
        mfd <- matrix.manifold('spd','LogEuclidean',c(d,d))
        
        for(n in c(1,10))
        {
            mu <- diag(rep(1,d)) #rmatrix.spd_LogCholesky(mfd,n=1,sig=0.1,drop=T)
            
            S <- rsym(d,n,drop=F)
            expect_true(all(dim(S)==c(mfd$dim,n)))
            
            S <- rmatrix.spd_LogEuclidean(mfd,n,drop=F,mu=mu)
            expect_true(all(sapply(1:n, function(i){
                is.spd(S[,,i])
            })))
            V <- rtvecor.spd_LogEuclidean(mfd,n=n,sig=0.1,drop=F)
            expect_true(all(sapply(1:n, function(i){
                is.sym(V[,,i])
            })))
            
            mu <- rmatrix.spd_LogEuclidean(mfd,n=1)
            V <- center.matrices(V)
            Q <- array(0,c(d,d,n))
            for(i in 1:n)
            {
                Q[,,i] <- rie.exp.spd_LogEuclidean(mfd,mu,V[,,i])
            }
            Q.mu <- frechet.mean.spd_LogEuclidean(mfd,Q)
            
            expect_true(frobenius.norm(mu-Q.mu) < eps)
        }
        
    }
    
    
})