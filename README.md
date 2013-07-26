
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

Enjoy! :)

