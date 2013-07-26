#####################################
## Generic Monte Carlo EM function ##
#####################################

#' Function for executing an instance of the Monte Carlo EM algorithm
#'
#' The function \code{\link{MCEM}} takes an \code{update} function, observed
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
#' When using \code{MCEM} it is expected that \code{update} will be a stochastic
#' function, and thus return (possibly) different values even if called with 
#' the same input arguments.
#' @param max.iter maximum number of iterations to run the algorithm
#' @param monitor whether to perform advanced convergence monitoring using
#' the methods in Baines, Xu and Wang (2013). If \code{FALSE} then basic
#' (but potentially unreliable) convergence monitoring using standard
#' relative or absolute tolerance is used. If \code{TRUE} then the estimate
#' is taken as the fixed point of the update mapping, which is approximated
#' via either a linear or spline fit to the EM path.  
#' @param ss.MC.reps the number of calls to \code{update} used to estimate
#' the Monte-Carlo variability in the update function.
#' @param mu.start starting value for mu
#' @param sigma.start starting value for sigma
#' @param iter.to.wait specifies the number of iterations to wait before checking
#' for convergence begins using a linear fit
#' @param cr the confidence level used in convergence checking
#' @param n.monitor.samples the number of monitoring samples to use
#' @param smooth whether to use a smooth function to approximate the update
#' mapping. If \code{TRUE} then a spline fit is used, otherwise a simple
#' linear fit is used.
#' @param penalty the penalty term used if a smooth function is used to approximate
#' the update mapping
#' @param B the bootstrap sample size used for approximating the update mapping
#' @param iter.to.wait.sp the number of iterations to wait before attempting
#' to compute a spline fit to the update mapping
#' @param spline.se.every the number of iterations at which the SE is recomputed
#' from the spline fit to the update mapping. Since this can be slow, it is not
#' recommended to compute this at each iteration.
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
#' @seealso \code{\link{EM}}
#' @examples
#' "MCEM.example.update" <- function(theta,y.obs,fixed,verbose)
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
#' ss.MC.reps <- 10
#' mu.start <- 0
#' sigma.start <- 1e5
#' iter.to.wait <- 5
#' cr <- 0.95
#' n.monitor.samples <- 10000
#' max.iter <- 2000 # 10 million
#' verbose <- FALSE # TRUE
#' print.every <- 100
#' smooth <- TRUE
#' penalty <- 1.0
#' B <- 100
#' spline.se.every <- 10000000000 # 100
#' 
#' t3 <- MCEM(y.obs=y.obs, fixed=fixed, theta.0=theta.0, update=EM.example.update,
#'            max.iter=max.iter, monitor=TRUE, ss.MC.reps=ss.MC.reps,
#'            mu.start=mu.start, sigma.start=sigma.start,
#'            smooth=smooth,B=B,penalty=penalty,spline.se.every=spline.se.every,
#'            iter.to.wait=iter.to.wait, cr=cr,n.monitor.samples=n.monitor.samples,
#'            print.every=print.every, tol=zero.tol, verbose=verbose)
#' print(t3)
#' @export
"MCEM" <- function(y.obs,theta.0,fixed,update,max.iter,monitor=TRUE,ss.MC.reps=10,
                   mu.start=0.0,sigma.start=1e5,iter.to.wait=5,cr=0.95,
                   n.monitor.samples=10000,smooth=FALSE,penalty=2,B=100,iter.to.wait.sp=50,
                   spline.se.every=100,
                   logLike=NULL,keep.paths=TRUE,append.paths=FALSE,tol=1e-10,
                   tol.type="relative",print.every=Inf,compute.DM=FALSE,verbose=FALSE)
{
  
  # Notes:
  # -- theta.0 (i.e., theta) can be a list
  # -- update should be a function, taking four arguments only: theta, y.obs, fixed and verbose
  # -- logLike should be a function taking two arguments only: theta and y.obs
  # -- smooth=TRUE means applying a smooth function to approximate G()
  # -- theta.bar and theta.bar.10 only appears when monitor=TRUE, not when smooth=TRUE
  
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
  
  ########
  # HACK -- Recode this when we decide exactly which summaries to keep for
  # monitoring MCEM convergence. 
  ########

  if (keep.paths){
    if (append.paths){
      # Keep paths via appending to current path:
      paths <- list("theta"=matrix(nrow=0,ncol=length.theta),"logLike"=NULL)
      if (monitor){
        # Modified slope-based estimator:
        paths$theta.tilde <- matrix(nrow=0,ncol=length.theta)
        # Long-run average estimator:
        paths$theta.bar <- matrix(nrow=0,ncol=length.theta)
        # Running mean (last-10) estimator:
        paths$theta.bar.10 <- matrix(nrow=0,ncol=length.theta)
        # SD(theta.tilde) 
        paths$tilde.sd <- matrix(nrow=0,ncol=length.theta)
        # QR(theta.tilde) 
        paths$tilde.qr <- matrix(nrow=0,ncol=length.theta)
        # SD(theta.bar) 
        paths$bar.sd <- matrix(nrow=0,ncol=length.theta)
        # QR(theta.bar) 
        paths$bar.qr <- matrix(nrow=0,ncol=length.theta)
        # SD(theta.bar.10) 
        paths$bar.10.sd <- matrix(nrow=0,ncol=length.theta)
        # QR(theta.bar.10) 
        paths$bar.10.qr <- matrix(nrow=0,ncol=length.theta)
      }
      if (smooth){
        # smooth spline estimator:
        paths$theta.spline <- matrix(nrow=0,ncol=length.theta)
        # SD(theta.spline)
        paths$spline.sd <- matrix(nrow=0,ncol=length.theta)
        # QR(theta.spline)
        paths$spline.qr <- matrix(nrow=0,ncol=length.theta)
      }
    } else {
      # Keep paths by making massive matrix and deleting unwanted rows at end:
      paths <- list("theta"=matrix(NA,nrow=max.iter,ncol=length.theta),"logLike"=NULL)
      rownames(paths$theta) <- ppaste("iter_",1:max.iter)
      if (is.function(logLike)){
        paths$logLike <- rep(NA,max.iter)
      }
      if (monitor){
        paths$theta.tilde  <- matrix(nrow=max.iter,ncol=length.theta)
        paths$theta.bar    <- matrix(nrow=max.iter,ncol=length.theta)
        paths$theta.bar.10 <- matrix(nrow=max.iter,ncol=length.theta)
        paths$tilde.sd  <- matrix(nrow=max.iter,ncol=length.theta)
        paths$tilde.qr  <- matrix(nrow=max.iter,ncol=length.theta)
        paths$bar.sd    <- matrix(nrow=max.iter,ncol=length.theta)
        paths$bar.qr    <- matrix(nrow=max.iter,ncol=length.theta)
        paths$bar.10.sd <- matrix(nrow=max.iter,ncol=length.theta)
        paths$bar.10.qr <- matrix(nrow=max.iter,ncol=length.theta)
      }
      if (smooth){
        paths$theta.spline <- matrix(nrow=max.iter,ncol=length.theta)
        paths$spline.sd <- matrix(nrow=max.iter,ncol=length.theta)
        paths$spline.qr <- matrix(nrow=max.iter,ncol=length.theta)
      }
    }
  }else {
    paths <- NULL
  }
  
  if (is.function(logLike)){
    ll.t <- logLike(theta.t,y.obs)
  }
  
  # EM-update approximation (for SEM and diagnostics):
  if (compute.DM && !keep.paths){
    stop("Computation of DM matrix requires 'keep.paths'=TRUE")
  }
  
  if (monitor){
    if (!is.list(mu.start)){
      mu.start.list <- vector("list",length.theta)
      for (i in 1:length.theta){
        mu.start.list[[i]] <- rep(mu.start[1+(i-1)%%length.theta],2) # TODO: Fix this...
      }
      mu.start <- mu.start.list
    }
    if (!is.list(sigma.start)){
      sigma.start.list <- vector("list",length.theta)
      for (i in 1:length.theta){
        sigma.start.list[[i]] <- diag(sigma.start[1+(i-1)%%length.theta],2)
      }
      sigma.start <- sigma.start.list
    }
    # Run pre-specified number of iterations to estimate MC-SD:
    if (verbose){
      cat("Estimating MC-SD...\n")
    }
    init <- mcem.monitor.initialize(theta.0=theta.0,y.obs=y.obs,fixed=fixed,EM.update=EM.example.update,iter.to.wait=iter.to.wait,ss.MC.reps=ss.MC.reps,verbose=verbose)    
    if (verbose){
      cat("Finished estimating MC-SD:\n")
      cat("MC-SD:\n")
      print(sqrt(init$ss.MC))
    }
    # Initialize MCEM convergence monitoring:
    if (verbose){
      cat("Initializing MCEM convergence monitoring...\n")
    }
    states <- mcem.monitor.start(theta.start=theta.0,mu.start=mu.start,sigma.start=sigma.start,ss.MC=init$ss.MC,cr=cr)
    if (verbose){
      cat("Completed initialization of MCEM convergence monitoring...\n")
    }
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
    
    if (monitor){
      # Perform MCEM convergence monitoring:
      monitor.t1 <- mcem.monitor.update(theta.t=theta.t,theta.t1=theta.t1,n.samples=n.monitor.samples,states=states)
      # Update state for future iterations:
      states <- monitor.t1$states
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
        if (monitor){
          # Need to compute these from the paths not just last state...
          monitor.t1$theta.bar <- apply(paths$theta,2,mean)
          iter.10 <- max(1,iter+1-9):(iter+1)
          b.hat <- extract.b.hat(monitor.t1,iter=iter+1)
          monitor.t1$theta.bar.10 <- apply(paths$theta[iter.10,,drop=FALSE],2,mean) # Note: Need theta.t1 to be appended first!
          bar.unc    <- bar.uncertainty(sigma.mc.sq=monitor.t1$states$ss.MC,b.hat=b.hat,n=iter+1,cr=cr)
          bar.10.unc <- bar.uncertainty(sigma.mc.sq=monitor.t1$states$ss.MC,b.hat=b.hat,n=length(iter.10),cr=cr)
          monitor.t1$bar.sd <- bar.unc$bar.sd
          monitor.t1$bar.qr <- bar.unc$bar.qr
          monitor.t1$bar.10.sd <- bar.10.unc$bar.sd
          monitor.t1$bar.10.qr <- bar.10.unc$bar.sd
          paths$theta.tilde <- append(paths$theta.tilde,monitor.t1$theta.tilde)
          paths$tilde.sd <- rbind(paths$tilde.sd,monitor.t1$tilde.sd)
          paths$tilde.qr <- rbind(paths$tilde.qr,monitor.t1$tilde.qr)          
          paths$theta.bar <- append(paths$theta.bar,monitor.t1$theta.bar)
          paths$theta.bar.10 <- append(paths$theta.bar.10,monitor.t1$theta.bar.10)
          paths$bar.sd <- rbind(paths$bar.sd,monitor.t1$bar.sd)
          paths$bar.qr <- rbind(paths$bar.qr,monitor.t1$bar.qr)          
          paths$bar.10.sd <- rbind(paths$bar.10.sd,monitor.t1$bar.10.sd)
          paths$bar.10.qr <- rbind(paths$bar.10.qr,monitor.t1$bar.10.qr)          
        }
        if (smooth){
          if (iter>=iter.to.wait.sp){
            spline.t1 <- mcem.spline.update(paths$theta,iter,penalty,cr,B,spline.se.every)
            paths$theta.spline <- append(paths$theta.spline,spline.t1$theta.spline)
            paths$spline.sd <- rbind(paths$spline.sd,spline.t1$spline.sd)
            paths$spline.qr <- rbind(paths$spline.qr,spline.t1$spline.qr)
          }     
        }
      }else {
        # Non-appending version:
        paths$theta[iter+1,] <- theta.t1
        if (is.function(logLike)){
          paths$logLike[iter+1] <- ll.t1
        }
        if (monitor){
          # Need to compute these from the paths not just last state...
          monitor.t1$theta.bar <- apply(paths$theta[1:(iter+1),,drop=FALSE],2,mean)
          if (verbose){
              cat("theta path:\n")
              print(paths$theta[1:(iter+1),])
          }
          iter.10 <- max(1,iter+1-9):(iter+1)
          b.hat <- extract.b.hat(monitor.t1,iter=iter+1)
          monitor.t1$theta.bar.10 <- apply(paths$theta[iter.10,,drop=FALSE],2,mean) # Note: Need theta.t1 to be appended first!
          bar.unc    <- bar.uncertainty(sigma.mc.sq=monitor.t1$states$ss.MC,b.hat=b.hat,n=iter+1,cr=cr)
          bar.10.unc <- bar.uncertainty(sigma.mc.sq=monitor.t1$states$ss.MC,b.hat=b.hat,n=length(iter.10),cr=cr)
          monitor.t1$bar.sd <- bar.unc$bar.sd
          monitor.t1$bar.qr <- bar.unc$bar.qr
          monitor.t1$bar.10.sd <- bar.10.unc$bar.sd
          monitor.t1$bar.10.qr <- bar.10.unc$bar.sd
          paths$theta.tilde[iter+1,] <- monitor.t1$theta.tilde
          paths$tilde.sd[iter+1,] <- monitor.t1$tilde.sd
          paths$tilde.qr[iter+1,] <- monitor.t1$tilde.qr
          paths$theta.bar[iter+1,] <- monitor.t1$theta.bar
          paths$theta.bar.10[iter+1,] <- monitor.t1$theta.bar.10
          paths$bar.sd[iter+1,] <- monitor.t1$bar.sd
          paths$bar.qr[iter+1,] <- monitor.t1$bar.qr
          paths$bar.10.sd[iter+1,] <- monitor.t1$bar.10.sd
          paths$bar.10.qr[iter+1,] <- monitor.t1$bar.10.qr
        }
        if (smooth){
          if (iter>=iter.to.wait.sp){
            spline.t1 <- mcem.spline.update(paths$theta,iter,penalty,cr,B,spline.se.every)
            paths$theta.spline[iter+1,] <- spline.t1$theta.spline
            paths$spline.sd[iter+1,] <- spline.t1$spline.sd
            paths$spline.qr[iter+1,] <- spline.t1$spline.qr
          }
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
    }else {
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
    if (monitor){
      # Slope-based estimate:
      paths$theta.tilde <- paths$theta.tilde[1:iter,,drop=FALSE]
      paths$tilde.sd <- paths$tilde.sd[1:iter,,drop=FALSE]
      paths$tilde.qr <- paths$tilde.qr[1:iter,,drop=FALSE]
      # Long-run mean estimate:
      paths$theta.bar <- paths$theta.bar[1:iter,,drop=FALSE]
      paths$bar.sd <- paths$bar.sd[1:iter,,drop=FALSE]
      paths$bar.qr <- paths$bar.qr[1:iter,,drop=FALSE]
      # Long-run mean estimate:
      paths$theta.bar.10 <- paths$theta.bar.10[1:iter,,drop=FALSE]
      paths$bar.10.sd <- paths$bar.10.sd[1:iter,,drop=FALSE]
      paths$bar.10.qr <- paths$bar.10.qr[1:iter,,drop=FALSE]
    }
    if (smooth){
      # Spline-based estimate:
      paths$theta.spline <- paths$theta.spline[1:iter,,drop=FALSE]
      paths$spline.sd <- paths$spline.sd[1:iter,,drop=FALSE]
      paths$spline.qr <- paths$spline.qr[1:iter,,drop=FALSE]
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
      
      
    
#####################################
##        Other subfunctions       ##
#####################################

"mcem.monitor.initialize" <- function(theta.0,y.obs,fixed,EM.update,iter.to.wait=5,ss.MC.reps=10,verbose=FALSE)
{
  # Run iter.to.wait iterations of EM...
  # Then repeat iter.to.wait-th step ss.MC.reps times, 
  # and use results to estimate sigma.MC:
  
  # Run EM without allowing convergence:
  initial.run <- EM(y.obs=y.obs,theta.0=theta.0,fixed=fixed,update=EM.update,max.iter=iter.to.wait,
                    logLike=NULL,keep.paths=FALSE,tol=0,tol.type="relative",print.every=Inf,compute.DM=FALSE,
                    verbose=verbose)
  
  # Requires theta is a vector:
  if (is.list(initial.run$theta)){
    stop("'theta' cannot be a list")
  }
  p <- length(as.numeric(initial.run$theta))
  theta.rep.list <- matrix(NA,nrow=ss.MC.reps,ncol=p)
  
  for (i in 1:ss.MC.reps){
    # Run 1 iteration:
    theta.rep.list[i,] <- as.numeric(EM(y.obs=y.obs,theta.0=initial.run$theta,fixed=fixed,update=EM.update,max.iter=1,
                                        logLike=NULL,keep.paths=FALSE,tol=0,tol.type="relative",print.every=Inf,compute.DM=FALSE,
                                        verbose=verbose)$theta)
  }
  
  ss.MC <- apply(theta.rep.list,2,var)
  
  return(list("ss.MC"=ss.MC))
}

"mcem.monitor.start" <- function(theta.start,mu.start,sigma.start,ss.MC,cr=0.95)
{
  # Set up the 'state' list:
  state <- list("mu.old"=mu.start,
                "sigma.old"=sigma.start,
                "ss.MC"=ss.MC,
                "p"=length(theta.start),
                "cr"=cr)
  
  return(state)
}

"mcem.monitor.update" <- function(theta.t,theta.t1,n.samples,states)
{
  # Dimensionality of theta?
  p <- states$p
  
  # New state (keep names):
  theta.tilde <- theta.t
  theta.tilde[1:p] <- NA
  tilde.sd <- tilde.qr <- rep(NA,p)
  cr <- states$cr ## ADDED 02/20/13
  
  for (i in 1:p){
    
    # Extract parameter-specific current states:
    sigma.tmp.old <- states$sigma.old[[i]]
    mu.tmp.old <- states$mu.old[[i]]
    ss.tmp.MC <- states$ss.MC[i]
    
    # Update linear approximation for alpha:
    Z.tmp <- theta.t[i]^(0:1)    
    sigma.tmp.new <- solve(solve(sigma.tmp.old)+Z.tmp%*%t(Z.tmp)/ss.tmp.MC[i])
    mu.tmp.new <- sigma.tmp.new%*%(solve(sigma.tmp.old)%*%mu.tmp.old+theta.t1[i]*Z.tmp/ss.tmp.MC[i])
    Sample.tmp <- try({rmvnorm(n=n.samples,mean=mu.tmp.new,sigma=sigma.tmp.new,method="chol")},silent=TRUE) # Cholesky for speed...
    if (class(Sample.tmp)=="try-error"){
      Sample.tmp <- rmvnorm(n=n.samples,mean=mu.tmp.new,sigma=sigma.tmp.new,method="svd") # Fall back...
    }
    # Compute \tilde estimates:
    tmp.theta.tilde <- Sample.tmp[,1]/(1-Sample.tmp[,2])
    tilde.sd[i] <- sd(tmp.theta.tilde)
    tilde.qr[i] <- diff(quantile(tmp.theta.tilde,probs=c(0.5*(1-cr),cr+0.5*(1-cr))))
    theta.tilde[i] <- median(tmp.theta.tilde)
    
    # Update the current states for parameter i:
    states$mu.old[[i]] <- mu.tmp.new
    states$sigma.old[[i]] <- sigma.tmp.new
    
  }
  
  return(list("theta.tilde"=theta.tilde,"tilde.sd"=tilde.sd,"tilde.qr"=tilde.qr,"states"=states))
}

"bar.uncertainty" <- function(sigma.mc.sq,b.hat,n,cr){
  # For last-10 version, just call with n=10...
  # Compute variance:
  bar.var <- 0.0
  for (j in 1:n){
    bar.var <- bar.var + sum(b.hat^(abs(j-c(1:n))))
  }
  bar.var <- sigma.mc.sq*bar.var/((1-(b.hat^2))*(n^2))
  # SD and normal-approximation QR:
  bar.sd <- sqrt(bar.var)
  bar.qr <- diff(qnorm(c(0.5*(1-cr),cr+0.5*(1-cr)),sd=bar.sd))
  return(list("bar.sd"=bar.sd,"bar.qr"=bar.qr))
}

"extract.b.hat" <- function(monitor,iter)
{
  p <- monitor$states$p
  bh <- rep(NA,p)
  for (i in 1:p){
    bh[i] <- monitor$states$mu.old[[i]][2]
  }
  return(bh)
}


"mcem.spline.update" <- function(paths.theta,iter,penalty,cr,B,spline.se.every){
  # function to use spline to fit the paths.theta
  p <- ncol(paths.theta)
  theta.spline <- rep(NA,p)
  spline.sd <- spline.qr <- rep(NA,p)
  
  for(k in 1:p){
    x <- paths.theta[1:(iter-1),k]
    y <- paths.theta[2:iter,k]
    ss <- smooth.spline(x=x,y=y,penalty=penalty)
    sol <- uniroot(f=G.root,interval=range(x),ss=ss)
    theta.spline[k] <- sol$root
    if (iter%%spline.se.every == 0){
      bout <-bootstrap(x,y,penalty,cr,B)
      spline.sd[k] <- bout$spline.sd
      spline.qr[k] <- bout$spline.qr
    }
  }
  return(list("theta.spline"=theta.spline,"spline.sd"=spline.sd,"spline.qr"=spline.qr))
}


"G.root" <- function(x,ss){
  return(predict(object=ss,x=x)$y-x)
}


"bootstrap" <- function(x,y,penalty,cr,B,print.every=10){
  # function to bootstrap to calculate the standard error and quantile range
  # B is the number of boostrap samples
  
  thetanew <- rep(0,B)
  nsize = length(x)
  i <- 0
  while(i<=B){
    resample <- sample(nsize,nsize,TRUE)
    xnew <- x[resample]
    ynew <- y[resample]
    ssnew <- try(smooth.spline(x=xnew,y=ynew,penalty=penalty),TRUE)
    if (inherits(ssnew,"try-error")) next
    solnew <- try(uniroot(f=G.root,interval=range(xnew),ss=ssnew),TRUE)
    if (inherits(solnew,"try-error")) next
    thetanew[i] <- solnew$root
    i <- i+1
    if (i%%print.every == 0){
      cat(paste("Finished bootsrap sample ",i,"/",B,"...\n",sep=""))
    }
  }
  spline.sd <- sd(thetanew)
  spline.qr <- diff(quantile(thetanew,probs=c(0.5*(1-cr),cr+0.5*(1-cr))))
  
  return(list("spline.sd"=spline.sd,"spline.qr"=spline.qr))
}
