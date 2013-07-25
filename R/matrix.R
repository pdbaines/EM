
#' @export
"spd.chol" <- function(X,lower=FALSE)
{
  C <- chol(X)
  if (lower){
    return(t(C))
  } 
  return(C)
}

#' @export
"svd.root" <- function(X,tol=1e-50,right=FALSE) # Works on Matrix objects too
{
  s <- svd(X)
  s$d[abs(s$d)<tol] <- tol 
  if (!right){
    ret <- s$u %*% diag(sqrt(s$d))
  } else {
    ret <- tcrossprod(diag(sqrt(s$d)), s$v)
  }
  return(ret)
}

#' @export
"svd.inverse.root" <- function(X,tol=1e-50) # Works on Matrix objects too
{
  s <- svd(X)
  s$d[abs(s$d)<tol] <- tol
  return(s$u %*% diag(1.0/sqrt(s$d)))
}

"svd.solve" <- function(X,tol=1e-50)
{
  s <- svd(X)
  s$d[abs(s$d)<tol] <- tol
  return(s$u %*% diag(1.0/s$d) %*% t(s$v))
}

# Courtesy of Alex Blocker (modified slightly):
# Functions to handle diagonal matrices much faster (array is very slow)
"diag_mat" <- function(x,n) {
    if (missing(n))
      n <- length(x)
    y <- matrix(0,n,n)
    y[1L + 0L:(n-1L) * (n+1L)] <- x
    return(y)
}

#' @export
"diag_ind" <- function(n) {
    1L + 0L:(n-1L)*(n+1L)
}

## MAKE GENERIC...
#' @export
"elt.locs" <- function(X){
  # Function such that:
  # X@x == X[elt.locs(X)]
  #
  if (length(grep(class(X),c("dsCMatrix","dgCMatrix")))>0){
    # Number of non-zero elts per column:
    nnz.pc <- diff(X@p)
    # Total number of non-zero elts:
    nnz.tot <- sum(nnz.pc)
    # Number of rows and columns:
    nr <- nrow(X)
    nc <- ncol(X)
    # Index locations of non-zero elts:
    ix <- rep(NA,nnz.tot)
    # Retrieve non-zero index locations:
    for (i in 1:nc){
      # Skip unless there are non-zero elts:
      if (nnz.pc[i]>0){
        is <- 1+X@p[i]            # Index start (switched to 1-indexing)
        ie <- is + nnz.pc[i] - 1  # Index end
        row.start <- (i-1)*nr     # Starting elt location for column i
        ix[is:ie] <- 1 + row.start + X@i[is:ie] # Switch to 1-indexing
      }
    }
    return(ix)
  } else
    stop("Invalid class '",class(X),"' in 'elt.locs'")
}

#' @export
"ix2im" <- function(ix,nr,nc)
{
  ## Index vector to index matrix
  ## i.e., (i-1)*nr + j |--> (i,j)
  ##   ==> k |--> (floor((k-1)/nr), k-floor((k-1)/nr))
  im <- floor((ix-1)/nr)
  ri <- ix - nr*im
  ci <- 1+im
  return(cbind(ri,ci))
}

## MAKE GENERIC...
#' @export
"matrix.flatten" <- function(X,diagonal=FALSE,triangular=FALSE,lower=TRUE,sparse=FALSE,diag=TRUE)
{
  ## Function to flatten a matrix into a vector.
  ## Returns a list with two components:
  ## -- elts : A vector of k elements (k<=n*m), for regular matrices k==n*m,
  ##             for Matrix class objects, k can be lower for Sparse matrices.
  ## -- locs : A vector of integers specifying the location of elts. i.e.,
  ##             X[locs] == elts.
  ##
  ## sparse is only listened to if diagonal and triangular are both FALSE...

  if (!is.logical(diagonal))
    stop("'diagonal' must be TRUE or FALSE")
  
  if (!is.logical(triangular))
    stop("'triangular' must be TRUE or FALSE")
  
  if (diagonal)
    return(list("elts"=diag(X),"locs"=diag_ind(ncol(X))))
  
  x.ind <- seq.int(1,length(X))

  if (!triangular)
    return(list("elts"=as.numeric(X),"locs"=x.ind))

  if (triangular){
    ## Triangular
    if (lower){
     locs <- x.ind[as.logical(lower.tri(X,diag=diag))]
    } else {
     locs <- x.ind[as.logical(upper.tri(X,diag=diag))]
    }
    return(list("elts"=X[locs],"locs"=locs))
  }

  if (sparse){
    if (is(X,"CsparseMatrix")){
      ## Only works for Csparse!
      elts <- X@x
      locs <- elt.locs(X)
      return(list("elts"=X@x,"locs"=locs))
    }
    stop("class '",class(X),"' not yet implemented in 'matrix.flatten'")
  }
}

#' @export
"tri.extract" <- function(X,lower=TRUE,diag=TRUE)
{
  if (lower){
   return(as.numeric(X[lower.tri(X,diag=diag)]))
  } else {
   return(as.numeric(X[upper.tri(X,diag=diag)]))
  }
}

#' @export
"matrix.expand" <- function(X,p,q,triangular=FALSE,lower=TRUE,diag=TRUE)
{
  if (!triangular){
    dim(X) <- c(p,q)
    return(X)
  }

  if (triangular){
    X <- tri.construct(X,lower=lower,diag=diag)
    return(X)
  }
}

#' @export
"tri.construct" <- function(x,lower=TRUE,diag=TRUE)
{
  # Check the size of the implied triangular matrix:
  N <- length(x)
  if (diag){
    q <- (sqrt(8*N+1)-1)/2
  } else {
    q <- (sqrt(8*N+1)+1)/2
  }
  # Must be an integer:
  if (abs(q-as.integer(round(q))) > .Machine$double.eps)
    stop(ppaste("Invalid length (",N,") in 'tri.construct'"))
  
  X <- matrix(0,q,q)
  if (lower){
    X[as.logical(lower.tri(X,diag=diag))] <- as.numeric(x) # Need to flatten to a vector for Matrix objects
  } else {
    X[as.logical(upper.tri(X,diag=diag))] <- as.numeric(x) # Need to flatten to a vector for Matrix objects
  }
  return(X)
}

#' @export
"is.square" <- function(X){
  if (nrow(X)==ncol(X)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @export
"is.symmetric" <- function(X,tol=.Machine$double.eps){
  if (sum(abs(X-t(X)))<tol){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @export
"is.spd" <- function(X,tol=.Machine$double.eps){
  if (!is.square(X)){
    warning("Non symmetric p.d. matrix: (non-square)")
    return(FALSE)
  }
  if (!is.symmetric(X,tol=tol)){
    warning("Non-symmetric p.d. matrix (non-symmetric)")
    return(FALSE)
  }
  if (!is.pd(X,tol=tol)){
    warning("Non-symmetric p.d. matrix (non-p.d.)")
    return(FALSE)
  }
  # All good:
  return(TRUE)
}

#' @export
"is.pd" <- function(X,tol=0.0)
{
  # Require for eigenvalues to be >tol, 
  # Select tol>0 for numerical stability... (no abs!)
  return(ifelse(all(eigen(X)$values>tol),TRUE,FALSE))
}

#' @export
"mirror.matrix" <- function(X,lowertoupper=TRUE,diag=FALSE)
{
  n <- nrow(X)
  if (ncol(X)!=n)
    stop("Non-square matrix in 'mirror.matrix'")

  if (n==1)
    return(X)

  if (diag){
    ilo <- 1
  } else {
    ilo <- 2
  }

  for (i in ilo:n){
    for (j in 1:ifelse(diag,i,i-1)){
      if (lowertoupper){
        X[j,i] <- X[i,j]
      } else {
        X[i,j] <- X[j,i]
      }
    }
  }

  return(X)
}

#' @export
"make.positive.definite" <- 
function (m, tol, check.method = "chol") 
{
    method = match.arg(check.method)
    if (!is.matrix(m)) 
        m = as.matrix(m)
    if (method == "eigen") {
        esv = eigen(m, only.values = TRUE)$values
        if (is.complex(esv)) {
            stop("Input matrix has complex eigenvalues!")
        }
        if (missing(tol)) 
            tol = max(dim(m)) * max(abs(esv)) * .Machine$double.eps
        if (sum(esv > tol) == length(esv)){
            return(m)
        } else {
            # Not p.d. patch up, recall with eigenvectors as well:
            es = eigen(m)
            esv = es$values
            delta = 2 * tol
            tau = pmax(0, delta - esv)
            dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
            return(m + dm)
        }
    }
    if (method == "chol") {
        val = try(chol(m), silent = TRUE)
        if (class(val) == "try-error"){
            # Not p.d. patch up, re-call with eigenvectors as well:
            es = eigen(m)
            esv = es$values
            d = max(dim(m))
            if (missing(tol)) 
              tol = d * max(abs(esv)) * .Machine$double.eps
            delta = 2 * tol
            tau = pmax(0, delta - esv)
            dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
            return(m + dm)
        } else {
           return(m)
        }
    }
}

