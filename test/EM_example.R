
# Start with clean slate:
rm(list=ls())

library(EM)

on.gauss <- FALSE

if (on.gauss){
  outdir <- "/home/pdbaines/Research/Convergence/"
} else {
  outdir <- "~/Dropbox/"  
}

"EM.example.update" <- function(theta,y.obs,fixed,verbose)
{
  y1 <- y.obs$y1
  y4 <- y.obs$y4
  n <- y.obs$n
  
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
fixed <- list()
theta.0 <- 0.5
max.iter <- 1000
verbose <- FALSE

cat("Running vanilla-EM...\n")
t1 <- EM(y.obs=y.obs, fixed=fixed, theta.0=theta.0, update=EM.example.update, max.iter=max.iter, verbose=verbose)
t1
cat("Done.\n")

cat("Making plots...\n")

plot(t1$paths$theta,ylab="theta^{(t)}",xlab="iteration",
main="EM Algorithm Path")
abline(h=t1$theta,col="red",lwd=1.5)

plot(diff(t1$paths$theta),ylab="theta^{(t+1)}-theta^{(t)}",xlab="iteration",
main="EM Algorithm Path: Increments")
abline(h=0,col="red",lwd=1.5)

plot(log(abs(diff(t1$paths$theta))),ylab="log(|theta^{(t+1)}-theta^{(t)}|)",xlab="iteration",
main="EM Algorithm Path: Log-Increments")

# Linear convergence:
updateplot(t1,main="EM Updates: theta^{(t+1)} vs. theta^{(t)}",addline=TRUE)

cat("done. :)\n")

