#########################
## Generic EM function ##
#########################

#' Function for executing an instance of the EM algorithm
#'
#' The function \code{\link{EM}} takes an \code{update} function, observed
#' data \code{y.obs}, starting values \code{theta.0}, any fixed quantities
#' \code{fixed} (see below) and other arguments and runs the algorithm
#' until a stopping rule is satisfied or a maximum number of iterations is 
#' reached.  
#' @param y.obs observed data (any format)
#' @param theta.0 initial value for parameter
#' @param fixed any additional fixed quantities (mainly for computational efficiency)
#' @param update function that outputs theta^{(t+1)} given theta^{(t)}. Must take
#' take arguments \code{theta.t}, \code{y.obs}, \code{fixed} and \code{verbose}
#' (in that order) i.e., it will be called as:
#'  \code{theta.t1 <- update(theta.t,y.obs,fixed,verbose)}
#' @param max.iter maximum number of iterations to run the algorithm
#' @param logLike (optional) function to compute the log-likelihood (or log-posterior)
#' at each iteration. Must take only two arguments: \code{theta.t} and 
#' \code{y.obs}. Will be called as e.g., \code{logLike(theta.t,y.obs)}
#' @param keep.paths whether to store the value of the parameter at each
#' iteration, or just return the final converged value. 
#' @param append.paths only used if \code{keep.paths=TRUE}. 
#' if \code{TRUE} then parameter values are appended after
#' each iteration. This is fast initially, but if lots of iterations are required
#' it will slow down. If \code{FALSE}, a full block of size \code{max.iter}
#' is allocated to store the parameter values, which can be unnecessary if the
#' algorithm converges quickly.
#' @param tol the tolerance used when checking for convergence at each iteration
#' @param tol.type the type of stopping rule used. Options are \code{relative}
#' and \code{absolute}.
#' @param print.every specifies the interval at which status updates are printed
#' to the screen e.g., \code{print.every=2} will print to screen after every
#' two iterations are completed.
#' @param compute.DM now defunct
#' @param verbose integer controling the level of verbosity of the function     
#' @details Lots of stuff
#' @return Lots of stuff
#' @examples
#' # From (ref)
#' "EM.example.update" <- function(theta,y.obs,fixed,verbose)
#' {
#'  y1 <- y.obs$y1
#'  y4 <- y.obs$y4
#'  n <- y.obs$n
#'  sigma.mc <- fixed$sigma.mc
#' 
#'  "mstep" <- function(y12, y11, y4, n)
#'  {
#'    return((y12+y4)/(n-y11))
#'  }
#'  "estep" <- function(psi_current, y1)
#'  {
#'    y11 <- (y1/2)/((1/2) + (psi_current/4))
#'    y12 <- y1 - y11
#'    return(list("y11"=y11,"y12"=y12))
#'  }
#' 
#'  E.vals <- estep(psi_current=theta, y1=y1)
#'  psi <- rnorm(n=1,mean=mstep(y12=E.vals$y12, y11=E.vals$y11, y4, n),sd=sigma.mc)
#' 
#'  return(psi)
#' }
#' 
#' y.obs <- list("y1"=125,"y2"=18,"y3"=20,"y4"=34,"n"=125+18+20+34)
#' fixed <- list("sigma.mc"=0)
#' theta.0 <- 0.5
#' max.iter <- 1000
#' verbose <- FALSE
#' t1 <- EM(y.obs=y.obs, fixed=fixed, theta.0=theta.0, update=EM.example.update, max.iter=max.iter, verbose=verbose)
#' t1
#' @export
"EM" <- function(y.obs,theta.0,fixed,update,max.iter,logLike=NULL,keep.paths=TRUE,append.paths=FALSE,tol=1e-10,tol.type="relative",print.every=Inf,compute.DM=FALSE,verbose=FALSE)
{

  # Notes:
  # -- theta.0 (i.e., theta) can be a list
  # -- update should be a function, taking fourarguments: theta, y.obs, fixed, verbose
  # -- logLike should be a function taking two arguments only: theta and y.obs


  if (verbose){
    cat("Observed data:\n\n")
    print(y.obs)
    cat("Fixed quantities:\n\n")
    print(fixed)
    cat("Starting state:\n\n")
    print(theta.0)
  }

  ret <- list("theta"=NULL,"iterations"=NULL,"logLike"=NULL)
  iter <- 0
  theta.t <- theta.0
  length.theta <- length(theta.t)


  # If the user wants to keep full paths, then two versions for storing paths
  # appending-- very slow but light on memory
  # pre-allocation -- fast but can use unecessary
  # TODO: implement smaller-block pre-alllocation + appending to adapt to size
  # without either extreme.

  if (keep.paths){
    if (append.paths){
      # Keep paths via appending to current path:
      paths <- list("theta"=matrix(nrow=0,ncol=length.theta),"logLike"=NULL)
    } else {
      # Keep paths by making massive matrix and deleting unwanted rows at end:
      paths <- list("theta"=matrix(NA,nrow=max.iter,ncol=length.theta),"logLike"=NULL)
      rownames(paths$theta) <- ppaste("iter_",1:max.iter)
      if (is.function(logLike)){
        paths$logLike <- rep(NA,max.iter)
      }
    }
  } else {
    paths <- NULL
  }

  # If available, compute log-likelihood of initial state:
  if (is.function(logLike)){
    ll.t <- logLike(theta.t,y.obs)
  }

  # EM-update approximation (for SEM and diagnostics):
  if (compute.DM && !keep.paths){
    stop("Computation of DM matrix requires 'keep.paths'=TRUE")
  }
  if (verbose)
    cat("Beginning EM algorithm...\n")

  converged <- FALSE

  while (iter < max.iter)
  {

    if (verbose){
      cat("theta^(t):\n")
      print(theta.t)
    }

    # Update:
    theta.t1 <- update(theta.t,y.obs,fixed,verbose)

    if (verbose){
      cat("theta^(t+1):\n")
      print(theta.t1)
    }

    if (is.function(logLike)){
      ll.t1 <- logLike(theta.t1,y.obs)
    }

    # Store the paths?
    if (keep.paths){
      if (append.paths){
        # Appending version:
        paths$theta   <- rbind(paths$theta,theta.t1)
        rownames(paths$theta)[nrow(paths$theta)] <- ppaste("iter_",iter)
        if (is.function(logLike)){
          paths$logLike <- append(paths$logLike,ll.t1)
        }
      } else {
        # Non-appending version:
        paths$theta[iter+1,] <- theta.t1
        if (is.function(logLike)){
          paths$logLike[iter+1] <- ll.t1
        }
      }
    }

    # Increment the iteration number
    iter <- iter + 1

    if (verbose)
      cat(ppaste("Finished iteration ",iter,"... checking for convergence...\n"))

    # Convergence check either on parameters, or, preferably, on the ll directly:
    if (is.function(logLike)){
      if (converged(ll.t1,ll.t,tol=tol,type=tol.type)){
        converged <- TRUE
        if (verbose)
          cat("\n === CONVERGED === \n\n")
        break
      }
    } else {
      if (converged(theta.t1,theta.t,tol=tol,type=tol.type)){
        converged <- TRUE
        if (verbose)
          cat("\n === CONVERGED === \n\n")
        break
      }
    }

    # Update:
    theta.t <- theta.t1
    if (is.function(logLike))
      ll.t <- ll.t1

    # Print to the user?
    pc <- (iter%%print.every==0)
    if (!is.na(pc) && pc){
      cat(ppaste("Finished iteration ",iter,"\n"))
    }

  }

  if (!converged)
    warning("Algorithm did not converge")

  if (!append.paths){
    # Need to trim extra rows:
    paths$theta <- paths$theta[1:iter,,drop=FALSE]
    if (is.function(logLike)){
      paths$logLike <- paths$logLike[1:iter]
    }
  }

  # Set up the return object:
  ret$converged  <- converged
  ret$iterations <- iter
  ret$paths      <- paths
  ret$theta      <- theta.t1
  if (is.function(logLike))
    ret$logLike <- ll.t1
  
  DM.path <- NULL
  if (compute.DM){
    ######
    # Optionally compute the DM matrix.
    # NOTE: Requires MLE, so must be computed post-hoc, and have keep.paths=TRUE
    cat("Computing DM matrix...\n")
    for (tmp.iter in 1:iter){
      DM.path[[tmp.iter]] <- try(DM.estimate(theta.mle=ret$theta,theta.t=ret$paths$theta[tmp.iter,],y.obs=y.obs,update=update,fixed=fixed,verbose=verbose),silent=TRUE)
      pc <- (tmp.iter%%print.every==0)
      if (!is.na(pc) && pc){
        cat(ppaste("Finished computing DM matrix for iteration ",tmp.iter,"\n"))
      }
    }
    cat("Finished computing DM matrix.\n")
    ######
  } else {
  }

  # Append DM information:
  ret$DM.path <- DM.path

  return(ret)
}

# Replaced with version in 'misc.R':
#
# "converged" <- function(theta.t1,theta.t,tol=1e-20,type="relative")
# {
#   if (type=="relative"){
# #Relative:
#     if (any(theta.t==0 | !is.finite(theta.t))){
# #Have to resort back to absolute here:
#       rats <- theta.t1
#     } else {
#       rats <- (theta.t1-theta.t)/theta.t
#     }
#   } else {
# #Absolute:
#     rats <- (theta.t1-theta.t)
#   } 
#   if (sqrt(sum(rats^2)) <= tol){
#     return(TRUE)
#   } else {
#     return(FALSE)
#   }
# }
