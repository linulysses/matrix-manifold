#' Transform a lower triangular matrix into a vector representation, diagonal by diagonal, starting from the main diagonal
#' @keywords internal
sym.to.vec.atomic <- function(S)
{
    if(is.vec(S) && length(S)==1) return(S)
    
    d <- dim(S)[1]
    v <- rep(0,d*(d+1)/2)
    v[1:d] <- diag(S)
    p <- d
    for(i in 2:d)
        for(j in 1:(d-i+1))
        {
            p <- p + 1
            v[p] <- S[i+j-1,j]
        }
    attr(v,'format') <- 'sym.in.vec'
    return(v)         
}

#' Transform a lower triangular matrix into a vector representation, diagonal by diagonal, starting from the main diagonal
#' @keywords internal
sym.to.vec <- function(S,drop=F)
{
    if(is.matrix(S)) #return(sym.to.vec.atomic(S))
        S <- array(S,dim=c(dim(S),1))
    
    if(is.array(S))
    {
        n <- count.sym.matrix(S)
        d <- dim(S)[1]
        V <- matrix(0,d*(d+1)/2,n)
        for(i in 1:n)
        {
            V[,i] <- sym.to.vec.atomic(S[,,i])
        }
        if(drop && n==1) V <- V[,1] 
        attr(V,'format') <- 'sym.in.vec'
        return(V)
    }
    else if(is.list(S))
    {
        lapply(S,sym.to.vec.atomic)
    }
}


#' Transform a vector created by \code{sym.to.vec.atomic} or \code{sym.to.vec} into a symmetric matric
#' @keywords internal
vec.to.sym.atomic <- function(v)
{
    d <- round(-0.5+sqrt(2*length(v)+0.25))
    S <- matrix(0,d,d)
    p <- 0
    for(i in 1:d)
        for(j in 1:(d-i+1))
        {
            p <- p + 1
            S[i+j-1,j] <- v[p]
            S[j,i+j-1] <- v[p]
        }
    attr(S,'format') <- 'sym.in.matrix'
    return(S)         
}

#' recover a (set of) lower triangular matrix from its vector representation
#' @keywords internal
vec.to.sym <- function(v,drop=F)
{
    if(is.vec(v))
    {
        v <- as.matrix(v)
        attr(v,'format') <- 'sym.in.vec'
    }
        
    
    if(is.matrix(v))
    {
        stopifnot(attr(v,'format') == 'sym.in.vec')
        
        n <- ncol(v)
        d <- round(-0.5+sqrt(2*length(v[,1])+0.25))
        S <- array(0,c(d,d,n))
        for(i in 1:n)
        {
            S[,,i] <- vec.to.sym.atomic(v[,i])
        }
        if(drop && n == 1) return(S[,,1])
        else return(S)
    }
    else if(is.list(v))
    {
        lapply(v,vec.to.sym.atomic)
    }
}

#' Check whether a symmetric matrix is represented by its vector form
#' @keywords internal
is.sym.in.vec.form <- function(S)
{
    if(is.list(S))
    {
        return(is.sym.in.vec.form(S[[1]]))
    }
    else if(is.array(S) || is.matrix(S))
    {
        return(attr(S,'format') == 'sym.in.vec')
    }
    else if(is.vec(S) && length(S) == 1) return(TRUE)
    else return(FALSE)
}

#' Check whether a matrix is represented by its vector form
#' @keywords internal
is.sym.in.matrix.form <- function(S)
{
    if(is.list(S))
    {
        return(is.sym.in.matrix.form(S[[1]]))
    }
    else if(!is.null(attr(S,'format')))
    {
        return(attr(S,'format') == 'sym.in.matrix')
    }
    else if(is.matrix(S))
    {
        return(is.sym(S))
    }
    else if(is.array(S))
    {
        return(is.sym.in.matrix.form(S[,,1]))
    }
    else return(FALSE)
}

#' transform a matrix into its matrix form
#' @keywords internal
sym.to.matrix.form <- function(S)
{
    if(is.sym.in.matrix.form(S)) return(S)
    else if(is.sym.in.vec.form(S))
    {
        return(vec.to.sym(S))
    }
    else stop('unrecoginzed format of S')
}

#' transform a matrix into its vector form
#' @keywords internal
sym.to.vec.form <- function(S)
{
    if(is.sym.in.vec.form(S)) return(S)
    else if(is.sym.in.matrix.form(S))
    {
        return(sym.to.vec(S))
    }
    else stop('unrecoginzed format of S')
}

#' get the matrix dimensions 
#' @keywords internal
get.sym.matrix.dim <- function(v)
{
    if(is.list(v)) return(get.sym.matrix.dim(v[[1]]))
    
    if(is.sym.in.vec.form(v))
    {
        if(is.matrix(v)) n <- nrow(v)
        else if(is.vec(v)) n <- length(v)
        else stop('unrecognized v')
        d <- round(-0.5+sqrt(2*n+0.25))
        return(d)
    }
    else if(is.sym.in.matrix.form(v))
    {
        if(is.array(v)) v <- v[,,1]
        return(nrow(v))
    }
    if(is.matrix(v)) return(nrow(v))
    else if(is.vec(v))
    {
        n <- length(v)
        d <- round(-0.5+sqrt(2*n+0.25))
        return(d)
    }
    else stop('unrecognized v')

}

#' count the number of matrices represented by an array or a list
#' @keywords internal
count.sym.matrix <- function(S)
{
    format <- attr(S,'format')
    if(!is.null(format) && format == 'sym.in.vec') return(ncol(S))
    else if(!is.null(format) && format == 'sym.in.matrix')
    {
        if(is.matrix(S)) return(1)
        else if(is.array(S)) return(dim(S)[3])
    }
    else if(is.array(S)) return(dim(S)[3])
    else if(is.matrix(S) && is.sym(S)) return(1)
    else if(is.list(S)) return(length(S))
    else stop('unrecognized format of S')
}
