context("check routines in symmetric.matrix.R")

test_that("check routines in symmetric.matrix.R", {
    
    eps <- 1e-12
    
    for(d in c(1,3))
    {
        for(n in c(1,10))
        {
            L <- make.sym(rmatrix(d,d))
            v <- sym.to.vec(L,drop=T)
            U <- vec.to.sym(v,drop=T)
            expect_true(all(L==U))
            
            S <- gen.sym(d,n,drop=F)
            V <- sym.to.vec(S,drop=F)
            expect_true(all(dim(V)==c(round(d*(d+1)/2),n)))
            R <- vec.to.sym(V,drop=F)
            expect_true(all(S==R))
            
            expect_true(is.sym.in.vec.form(V))
            expect_false(is.sym.in.vec.form(S))
            
            expect_true(is.sym.in.matrix.form(S))
            expect_false(is.sym.in.matrix.form(V))
            
            S <- make.sym(rmatrix(d,d))
            expect_true(is.sym.in.matrix.form(S))
            
            S1 <- sym.to.matrix.form(S)
            expect_true(all(S==S1))
            
            S <- gen.sym(d,n,drop=F)
            S1 <- sym.to.matrix.form(S)
            expect_true(all(S==S1))
            
            V <- sym.to.vec.form(S)
            S1 <- sym.to.matrix.form(V)
            expect_true(all(S==S1))
            
            m <- count.sym.matrix(S)
            expect_true(m==n)
        }
        
        
    }
    
    
})