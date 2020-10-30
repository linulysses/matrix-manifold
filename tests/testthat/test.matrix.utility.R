context("check routines in matrix.utility.R")

test_that("check routines in matrix.utility.R", {
    
    eps <- 1e-12
    
    for(d in c(1,3))
    {
        for(n in c(1,10))
        {
            L <- lower.part(rmat(d,d),strict=F)
            v <- lower.to.vec(L,drop=T)
            U <- vec.to.lower(v,drop=T)
            expect_true(all(L==U))
            
            S <- rltmatrix(d,n,drop=F)
            V <- lower.to.vec(S,drop=F)
            expect_true(all(dim(V)==c(round(d*(d+1)/2),n)))
            R <- vec.to.lower(V,drop=F)
            expect_true(all(S==R))
        }
    }
    
    for(d in c(1,3))
    {
        
        # check is.sym
        A <- rmat(d,d)
        S <- A + t(A)
        expect_true(is.sym(S))
        if(d > 1)  expect_false(is.sym(A))
        
        # check make.sym
        expect_true(is.sym(make.sym(A)))
        
        # check is.spd
        P <- make.sym(S %*% S)
        expect_true(is.spd(P))
        
        NP <- P
        NP[1,1] <- -1
        expect_false(is.spd(NP))
        
        # check is.vec
        expect_true(is.vec(runif(d)))
        expect_true(is.vec(rmat(1,d)))
        expect_true(is.vec(rmat(d,1)))
        if(d > 1) expect_false(is.vec(rmat(d,d)))
        
        # check is.lower
        L <- rmat(d,d)
        if(d > 1) expect_false(is.lower(L))
        for(i in 1:d)
            for(j in 1:d)
                if(j > i) L[i,j] <- 0
        expect_true(is.lower(L))
        
        L <- lower.part(rmat(d,d),strict=T)
        expect_true(is.lower(L))
        expect_true(all(diag(L)==0))
        L <- lower.part(rmat(d,d),strict=F)
        expect_false(all(diag(L)==0))
        
        U <- upper.part(rmat(d,d),strict=T)
        expect_true(is.lower(t(U)))
        U <- upper.part(rmat(d,d),strict=F)
        expect_false(all(diag(U)==0))
        
        # check frobenius
        A <- rmat(d,d)-0.5
        B <- rmat(d,d)-0.5
        
        expect_true(abs(sum(A*B)-frobenius(A,B)) < eps)
        expect_true(abs(sum(A*A)-frobenius.norm(A)^2) < eps)
    }
    
    # check half
    A <- rmat(d,d)
    H <- half(A)
    expect_true(all(strict.lower.part(t(H))==0))
    expect_true(all(diag(H) == diag(A)/2))
    
    # check operations on diagonal matrices
    D <- diag(diag(rmat(d,d)))
    ED <- expm.diag(D)
    expect_true(all(diag(ED)== exp(diag(D))))
    diag(ED) <- 0
    expect_true(all(ED==0))
    expect_true(frobenius.norm(logm.diag(D)-expm::logm(D)) < eps)
    E <- diag(diag(rmat(d,d)))
    expect_true(all(mul.diag(D,E)==D %*% E))
    expect_true(all(inv.diag(D)==solve(D)))
    
    # trace
    A <- rmat(d,d)
    expect_true(trace(A) == matrixcalc::matrix.trace(A))
    
})