
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
  t2 <- EM(y.obs=y.obs, fixed=fixed, theta.0=theta.0, update=EM.example.update,
           max.iter=max.iter, print.every=print.every, tol=zero.tol, verbose=verbose)
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
plot(rs,type="l",ylim=t1$theta+c(-0.0001,0.0001),ylab="theta",
     main="Monte Carlo EM: Long-run Average",
     xlab="iteration",sub=paste("sigma.mc=",sigma.mc,sep=""))
abline(h=mean(t2$paths$theta),col="red")
abline(h=mean(t1$theta),col="blue")
legend("topright",legend=c("MLE","Long-run Average"),col=c("blue","red"),lty=1)
dev.off()

pdf(paste(outdir,"MCEM_example_mapping.pdf",sep=""))
plot(x=t2$paths$theta[1:(t2$iterations-1)],
     y=t2$paths$theta[2:t2$iterations],
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
max.iter <- 2000 # 10 million
verbose <- FALSE # TRUE
print.every <- 100
smooth <- TRUE
penalty <- 1.0
B <- 100
spline.se.every <- 10000000000 # 100

cat(paste("Running MCEM with tilde and spline estimates for ",max.iter," iterations...\n",sep=""))
t3 <- MCEM(y.obs=y.obs, fixed=fixed, theta.0=theta.0, update=EM.example.update, 
           max.iter=max.iter, monitor=TRUE, ss.MC.reps=ss.MC.reps, 
           mu.start=mu.start, sigma.start=sigma.start,
           smooth=smooth,B=B,penalty=penalty,spline.se.every=spline.se.every,
           iter.to.wait=iter.to.wait, cr=cr,n.monitor.samples=n.monitor.samples,
           print.every=print.every, tol=zero.tol, verbose=verbose)
cat("done.\n")
          
t3$paths$theta
t3$paths$theta.tilde
t3$paths$tilde.sd
t3$paths$tilde.qr

cat("Making plots...\n")

t3$theta.tilde <- t3$paths$theta.tilde[length(t3$paths$theta.tilde)]
t3$theta.bar <- t3$paths$theta.bar[length(t3$paths$theta.bar)]
t3$theta.bar.10 <- t3$paths$theta.bar.10[length(t3$paths$theta.bar.10)]
t3$theta.spline <- t3$paths$theta.spline[length(t3$paths$theta.spline)]

ests <- c("MLE"=t1$theta, # true MLE
         "theta.tilde"=t3$theta.tilde ,
         "theta.bar"=t3$theta.bar,
         #"theta.bar.10"=t3$theta.bar.10,
         "theta.spline"=t3$theta.spline)

print(ests)

#B <- 10
theta.spline <- t3$paths$theta.spline
iter <- max.iter
k <- 1
x <- t3$paths$theta[1:(iter-1),k]
y <- t3$paths$theta[2:iter,k]
thetanew <- rep(0,B)
nsize <- length(x)
i <- 0
ss <- vector("list",B)
sols <- vector("list",B)
theta.s.vec <- rep(NA,B)
data.spline <- smooth.spline(x=x,y=y,penalty=penalty)
print.progress <- FALSE
while (i<=B){
  if (print.progress){
    cat("Resampling data...\n")
  }
  resample <- sample(nsize,nsize,TRUE)
  xnew <- x[resample]
  ynew <- y[resample]
  if (print.progress){
    cat("Fitting spline...\n")
  }
  ss[[i+1]] <- try(smooth.spline(x=xnew,y=ynew,penalty=penalty),TRUE)
  if (inherits(ss[[i+1]],"try-error")){
    cat("Error in spline fit!\n")
    warning("Error in spline fit")
    next
  } 
  sols[[i+1]] <- try(uniroot(f=G.root,interval=range(xnew),ss=ss[[i+1]]),TRUE)
  if (inherits(sols[[i+1]],"try-error")){
      cat("Error finding root!\n")
      warning("Error finding root")
      next
  }
  theta.s.vec[i] <- sols[[i+1]]$root
  i <- i+1
  if (i%%print.every == 0){
    cat(paste("Finished bootsrap sample ",i,"/",B,"...\n",sep=""))
  }
}

library(MASS)
truehist(theta.s.vec)
ss
sols

lims <- NULL # c(0.625,0.630)
plot(x=x,y=y,xlim=lims,ylim=lims)
for (i in 1:B){
  lines(ss[[i]],col="red")
}
lines(data.spline,col="green",lwd=2.0)
abline(a=0,b=1,col="blue",lwd=2.0)
abline(h=t1$theta,col="purple",lwd=2.0)
abline(v=t1$theta,col="purple",lwd=2.0)
points(x=t1$theta,y=t1$theta,pch=16,cex=1.5,col="purple")

matplot(cbind(t3$paths$theta,t3$paths$theta.tilde),type="l")

colvec <- rainbow(length(ests))
pdf(paste(outdir,"mcem_comparison.pdf",sep=""))
rs2 <- cbind(t3$paths$theta.bar,
             t3$paths$theta.tilde,
             #t3$paths$theta.bar.10,
             t3$paths$theta.spline)
matplot(rs2,type="l",ylab="theta",main="Monte Carlo EM: Comparison",
        ylim=range(rs2[floor(nrow(rs2)/4):nrow(rs2),],na.rm=TRUE),
        xlab="iteration",sub=paste("sigma.mc=",sigma.mc,sep=""),
        col=colvec[2:length(colvec)],lty=1)
abline(h=ests,col=colvec)
legend("topright",legend=names(ests),col=colvec,lty=1)
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

