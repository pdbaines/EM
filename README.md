
EM: An R Package
================

Author: Paul D. Baines
Original Date: 01/16/2011
Version Date: 01/30/2013

General tools for the EM Algorithm, including MCEM convergence monitoring.

Also see: http://pdbaines.github.io/EM/

To install:

R CMD INSTALL EM

To use:

(1) Write an EM-update function of the form:

  my.em.update.func <- function(theta.t,y.obs,fixed,verbose){
    # Do cool stuff...
    theta.t1 <- ... # Find update
    return(theta.t1)
  }

(2) Use the "EM" function to perform the EM algorithm with your update function,
    this includes all built-in convergence monitoring, options to save sample
    paths, store log-likehood paths (if the log-likelihood function is supplied), 
    and other useful stuff.

(2-alt) Use the "MCEM" function if your EM algorithm is a Monte Carlo EM algorithm.
     The MCEM function includes improved parameter estimation and convergence
     monitoring tools specifically tailored to MCEM algorithms.

Basic example:

"EM.example.update" <- function(theta,y.obs,fixed,verbose)
{
  y1 <- y.obs$y1
  y4 <- y.obs$y4
  n <- y.obs$n
  sigma.mc <- fixed$sigma.mc
 
  "mstep" <- function(y12, y11, y4, n)
  {
    return((y12+y4)/(n-y11))
  }
  "estep" <- function(psi_current, y1)
  {
    y11 <- (y1/2)/((1/2) + (psi_current/4))
    y12 <- y1 - y11
    return(list("y11"=y11,"y12"=y12))
  }
 
  E.vals <- estep(psi_current=theta, y1=y1)
  psi <- mstep(y12=E.vals$y12, y11=E.vals$y11, y4, n)
 
  return(psi)
}

y.obs <- list("y1"=125,"y2"=18,"y3"=20,"y4"=34,"n"=125+18+20+34)
fixed <- list("sigma.mc"=0)
theta.0 <- 0.5
t1 <- EM(y.obs=y.obs, fixed=fixed, theta.0=theta.0, update=EM.example.update, max.iter=100)
t1

Enjoy! :)

