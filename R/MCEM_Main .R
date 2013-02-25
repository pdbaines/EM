
#####################################
## Generic Monte Carlo EM function ##
#####################################

"MCEM" <- function(y.obs,theta.0,fixed,update,max.iter,monitor=TRUE,ss.MC.reps=10,
                   mu.start=0.0,sigma.start=1e5,iter.to.wait=5,cr=0.95,
                   n.monitor.samples=10000,logLike=NULL,
                   keep.paths=TRUE,append.paths=FALSE,tol=1e-10,tol.type="relative",
                   print.every=Inf,compute.DM=FALSE,verbose=FALSE)
{

  # Notes:
  # -- theta.0 (i.e., theta) can be a list
  # -- update should be a function, taking two arguments only: theta and y.obs
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
  
  if (keep.paths){
    if (append.paths){
      # Keep paths via appending to current path:
      paths <- list("theta"=matrix(nrow=0,ncol=length.theta),"logLike"=NULL)
      if (monitor){
        # Modified slope-based estimator:
        paths$theta.tilde <- matrix(nrow=0,ncol=length.theta)
        # Long-run estimator:
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
    }
  } else {
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
          paths$theta.bar <- append(paths$theta.tilde,monitor.t1$theta.bar)
          paths$theta.bar.10 <- append(paths$theta.tilde,monitor.t1$theta.bar.10)
          paths$bar.sd <- rbind(paths$bar.sd,monitor.t1$bar.sd)
          paths$bar.qr <- rbind(paths$bar.qr,monitor.t1$bar.qr)          
          paths$bar.10.sd <- rbind(paths$bar.10.sd,monitor.t1$bar.10.sd)
          paths$bar.10.qr <- rbind(paths$bar.10.qr,monitor.t1$bar.10.qr)          
        }
      } else {
        # Non-appending version:
        paths$theta[iter+1,] <- theta.t1
        if (is.function(logLike)){
          paths$logLike[iter+1] <- ll.t1
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
          paths$theta.tilde[iter+1,] <- monitor.t1$theta.tilde
          paths$tilde.sd[iter+1,] <- monitor.t1$tilde.sd
          paths$tilde.qr[iter+1,] <- monitor.t1$tilde.qr
          paths$theta.bar[iter+1,] <- monitor.t1$theta.bar
          paths$bar.sd[iter+1,] <- monitor.t1$bar.sd
          paths$bar.qr[iter+1,] <- monitor.t1$bar.qr
          paths$bar.10.sd[iter+1,] <- monitor.t1$bar.10.sd
          paths$bar.10.qr[iter+1,] <- monitor.t1$bar.10.qr
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
    bh[i] <- monitor$states$mu.old[[iter]][2]
  }
  return(bh)
}

