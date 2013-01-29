
"ppaste" <- function(...)
{
  paste(...,sep="")
}


"ldiff" <- function(x,y){
  if (length(x) != length(y)){
    stop("'x' and 'y' must have the same length")
  }
  if (is.list(x)){
    ret <- vector("list",length(x))
    names(ret) <- names(x)
    if (any(names(x)!= names(y))){
      warning("'x' and 'y' have different names")
    }
    for (i in 1:length(x)){
      ret[[i]] <- ldiff(x[[i]],y[[i]])
    }
  } else {
      ret <- (x-y)
  }
  return(ret)
}
      
"ldiff" <- function(x,y){
  if (length(x) != length(y)){
    stop("'x' and 'y' must have the same length")
  }
  if (is.list(x)){
    ret <- vector("list",length(x))
    names(ret) <- names(x)
    if (any(names(x)!= names(y))){
      warning("'x' and 'y' have different names")
    }
    for (i in 1:length(x)){
      ret[[i]] <- ldiff(x[[i]],y[[i]])
    }
  } else {
      ret <- (x-y)
  }
  return(ret)
}

"converged" <- function(theta.t1,theta.t,tol=1e-10,type="absolute")
{
  if (type=="relative"){
    # Relative:
    if (any(theta.t==0 | !is.finite(theta.t))){
      # Have to resort back to absolute here:
      rats <- theta.t1
    } else {
      rats <- (theta.t1-theta.t)/theta.t
    }
  } else {
    # Absolute:
    rats <- (theta.t1-theta.t)
  }
  if (all(abs(rats)<=tol)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

"NewtonRaphson" <- function(func,deriv.func,init,tol=1e-10,maxit=1000,simplify=TRUE,...)
{
  # Function to find x s.t. func(x)=0, where deriv.func(x)=dfunc(x)/dx
  converged.yet <- FALSE

  # Initial values:
  x.t <- init
  f.t <- func(x.t,...)

  for (i in 1:maxit){
    x.t <- x.t - solve(deriv.func(x.t,...))%*%f.t
    f.t <- func(x.t,...)
    if (converged(f.t,0.0,tol=tol,type="absolute")){
      converged.yet <- TRUE
      break
    }
  }

  if (!converged.yet)
    warning("Newton-Raphson algorithm failed to converge")

  if (simplify){
    return(x.t)
  } else {
    return(list("val"=x.t,"iterations"=maxit,"convergence"=converged.yet))
  }

}

