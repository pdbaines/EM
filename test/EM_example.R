
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
  psi <- rnorm(n=1,mean=mstep(y12=E.vals$y12, y11=E.vals$y11, y4, n),sd=sigma.mc)
  
  return(psi)
}

y.obs <- list("y1"=125,"y2"=18,"y3"=20,"y4"=34,"n"=125+18+20+34)
fixed <- list("sigma.mc"=0)
theta.0 <- 0.5
max.iter <- 1000
verbose <- FALSE

cat("Running vanilla-EM (no Monte Carlo Error)...\n")
t1 <- EM(y.obs=y.obs, fixed=fixed, theta.0=theta.0, update=EM.example.update, max.iter=max.iter, verbose=verbose)
t1
cat("Done.\n")

set.seed(5120417)

max.iter <- 10000 # 10 million
print.every <- 1000
zero.tol <- 0.0
sigma.mc <- 0.01
fixed <- list("sigma.mc"=sigma.mc)
# Never allow convergence:
cat(paste("Running EM with MC-error (no tilde estimates) for ",max.iter," iterations...\n",sep=""))
t2.time <- system.time({
  t2 <- EM(y.obs=y.obs, fixed=fixed, theta.0=theta.0, update=EM.example.update, max.iter=max.iter, print.every=print.every, tol=zero.tol, verbose=verbose)
})
cat("done. Timing:\n")
print(t2.time)

t2$theta
mean(t2$paths$theta)
t1$theta

cat("Making plots...\n")
rs <- cumsum(t2$paths$theta)/seq(along=t2$paths$theta) 
if (max.iter>1000000){
  pdf(paste(outdir,"mcem_long_run.pdf",sep=""))
} else {
  pdf(paste(outdir,"mcem_short_run.pdf",sep=""))
}
plot(rs,type="l",ylim=t1$theta+c(-0.0001,0.0001),ylab="theta",main="Monte Carlo EM: Long-run Average",
     xlab="iteration",sub=paste("sigma.mc=",sigma.mc,sep=""))
abline(h=mean(t2$paths$theta),col="red")
abline(h=mean(t1$theta),col="blue")
legend("topright",legend=c("MLE","Long-run Average"),col=c("blue","red"),lty=1)
dev.off()

pdf(paste(outdir,"MCEM_example_mapping.pdf",sep=""))
plot(x=t2$paths$theta[1:(t2$iterations-1)],
     y=t2$paths$theta[2:t2$iterations],cex.pt=0.6,
     xlab="theta^{(t)}",ylab="theta^{(t+1)}",
     main="MCEM Stochastic Update Mapping")
dev.off()

save.image(paste(outdir,"large_sim.RData",sep=""))

library(coda)

pdf(paste(outdir,"t2_path.pdf",sep=""))
plot(mcmc(t2$paths$theta))
dev.off()

ss.MC.reps <- 10
mu.start <- 0
sigma.start <- 1e5
iter.to.wait <- 5
cr <- 0.95
n.monitor.samples <- 10000
max.iter <- 1000 # 10 million
verbose <- FALSE # TRUE
print.every <- 100

cat(paste("Running MCEM with tilde estimates for ",max.iter," iterations...\n",sep=""))
t3 <- MCEM(y.obs=y.obs, fixed=fixed, theta.0=theta.0, update=EM.example.update, 
           max.iter=max.iter, monitor=TRUE, ss.MC.reps=ss.MC.reps, 
           mu.start=mu.start, sigma.start=sigma.start,
           iter.to.wait=iter.to.wait, cr=cr,n.monitor.samples=n.monitor.samples,
           print.every=print.every, tol=zero.tol, verbose=verbose)
cat("done.\n")
          
t3$paths$theta
t3$paths$theta.tilde
t3$paths$tilde.sd
t3$paths$tilde.qr

cat("Making plots...\n")

c("true.MLE"=t1$theta, # true MLE
  "Long.run.mean"=mean(t2$paths$theta), # long-run mean
   "Our.estimate"=t3$paths$theta.tilde[length(t3$paths$theta.tilde)]) # our estimate

matplot(cbind(t3$paths$theta,t3$paths$theta.tilde),type="l")

pdf(paste(outdir,"mcem_comparison.pdf",sep=""))
rs2 <- cbind(rs[1:t3$iterations],t3$paths$theta.tilde)
matplot(rs2,type="l",ylab="theta",main="Monte Carlo EM: Comparison",
        ylim=t1$theta+c(-0.0003,0.0003),
        xlab="iteration",sub=paste("sigma.mc=",sigma.mc,sep=""),
        col=c("red","orange"),lty=1)
abline(h=mean(t2$paths$theta[1:t3$iterations]),col="red")
abline(h=mean(t1$theta),col="blue")
abline(h=t3$paths$theta.tilde[length(t3$paths$theta.tilde)],col="orange")
legend("bottomright",legend=c("MLE","Long-run Average","tilde-estimate"),
       col=c("blue","red","orange"),lty=1)
dev.off()

plot(log(t3$paths$tilde.sd),type="l",main="MCEM Estimate SE",ylab="log(SD)",xlab="Iteration")
lines(log(t3$paths$bar.sd),col="red")
lines(log(t3$paths$bar.10.sd),col="blue")

plot(log(t3$paths$tilde.qr),type="l",main="MCEM Quantile Range",ylab="log(QR)",xlab="Iteration")
lines(log(t3$paths$bar.qr),col="red")
lines(log(t3$paths$bar.10.qr),col="blue")

cat("done. Saving image...\n")

# This will be large... :)
save.image(paste(outdir,"full_run_mcem.RData",sep=""))

cat("done. :)\n")

