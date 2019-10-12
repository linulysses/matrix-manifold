context("check routines in triangular.matrix.R")

test_that("check routines in triangular.matrix.R", {
    
    eps <- 1e-12
    
    for(d in c(1,3))
    {
        for(n in c(1,10))
        {
            L <- lower.part(rmatrix(d,d),strict=F)
            v <- lower.to.vec(L,drop=T)
            U <- vec.to.lower(v,drop=T)
            expect_true(all(L==U))
            
            S <- gen.lower(d,n,drop=F)
            V <- lower.to.vec(S,drop=F)
            expect_true(all(dim(V)==c(round(d*(d+1)/2),n)))
            R <- vec.to.lower(V,drop=F)
            expect_true(all(S==R))
        }
        
        
    }
    
    
})