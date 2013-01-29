
##
## Testsuite for rtnorm, dtnorm, qtnorm, ptnorm:
##

library(IEM)

n <- 10000
l <- 1.5
u <- 1.8
m <- 2.0
s <- 4.0

x1 <- rtnorm(n=n)
x2 <- rtnorm(n,l,u)
x3 <- rtnorm(n,lower=l,mean=m)
x4 <- rtnorm(n,l,u,m,s)

d1 <- dtnorm(x1)
d2 <- dtnorm(x2,l,u)
d3 <- dtnorm(x3,lower=l,mean=m)
d4 <- dtnorm(x4,l,u,m,s)

p1 <- ptnorm(x1)
p2 <- ptnorm(x2,l,u)
p3 <- ptnorm(x3,lower=l,mean=m)
p4 <- ptnorm(x4,l,u,m,s)

q1 <- qtnorm(p1)
q2 <- qtnorm(p2,l,u)
q3 <- qtnorm(p3,lower=l,mean=m)
q4 <- qtnorm(p4,l,u,m,s)


par(mfrow=c(2,2))
cex <- 0.8
cex.main <- 0.8
plot(x1,ylab="x",xlab="Sample id",main=ppaste("Sample from rtnorm (l=",-Inf,", u=",Inf,", m=",0,", sd=",1,")"),cex=cex,cex.main=cex.main)
plot(x2,ylab="x",xlab="Sample id",main=ppaste("Sample from rtnorm (l=",l,", u=",u,", m=",0,", sd=",1,")"),cex=cex,cex.main=cex.main) ; abline(h=c(l,u),col="red",lwd=1.5)
plot(x3,ylab="x",xlab="Sample id",main=ppaste("Sample from rtnorm (l=",l,", u=",Inf,", m=",m,", sd=",1,")"),cex=cex,cex.main=cex.main) ; abline(h=l,col="red",lwd=1.5)
plot(x4,ylab="x",xlab="Sample id",main=ppaste("Sample from rtnorm (l=",l,", u=",u,", m=",m,", sd=",s,")"),cex=cex,cex.main=cex.main) ; abline(h=c(l,u),col="red",lwd=1.5)

plot(y=d1,x=x1,xlab="x",ylab="f(x)",main=ppaste("Sample from rtnorm (l=",-Inf,", u=",Inf,", m=",0,", sd=",1,")"),cex=cex,cex.main=cex.main)
plot(y=d2,x=x2,xlab="x",ylab="f(x)",main=ppaste("Sample from rtnorm (l=",l,", u=",u,", m=",0,", sd=",1,")"),cex=cex,cex.main=cex.main)
plot(y=d3,x=x3,xlab="x",ylab="f(x)",main=ppaste("Sample from rtnorm (l=",l,", u=",Inf,", m=",m,", sd=",1,")"),cex=cex,cex.main=cex.main)
plot(y=d4,x=x4,xlab="x",ylab="f(x)",main=ppaste("Sample from rtnorm (l=",l,", u=",u,", m=",m,", sd=",s,")"),cex=cex,cex.main=cex.main)

plot(y=p1,x=x1,xlab="x",ylab="F(x)",main=ppaste("Sample from rtnorm (l=",-Inf,", u=",Inf,", m=",0,", sd=",1,")"),cex=cex,cex.main=cex.main)
plot(y=p2,x=x2,xlab="x",ylab="F(x)",main=ppaste("Sample from rtnorm (l=",l,", u=",u,", m=",0,", sd=",1,")"),cex=cex,cex.main=cex.main)
plot(y=p3,x=x3,xlab="x",ylab="F(x)",main=ppaste("Sample from rtnorm (l=",l,", u=",Inf,", m=",m,", sd=",1,")"),cex=cex,cex.main=cex.main)
plot(y=p4,x=x4,xlab="x",ylab="F(x)",main=ppaste("Sample from rtnorm (l=",l,", u=",u,", m=",m,", sd=",s,")"),cex=cex,cex.main=cex.main)

plot(y=q1,x=x1,xlab="x",ylab="F^{-1}F(x)",main=ppaste("Sample from rtnorm (l=",-Inf,", u=",Inf,", m=",0,", sd=",1,")"),cex=cex,cex.main=cex.main)
plot(y=q2,x=x2,xlab="x",ylab="F^{-1}F(x)",main=ppaste("Sample from rtnorm (l=",l,", u=",u,", m=",0,", sd=",1,")"),cex=cex,cex.main=cex.main)
plot(y=q3,x=x3,xlab="x",ylab="F^{-1}F(x)",main=ppaste("Sample from rtnorm (l=",l,", u=",Inf,", m=",m,", sd=",1,")"),cex=cex,cex.main=cex.main)
plot(y=q4,x=x4,xlab="x",ylab="F^{-1}F(x)",main=ppaste("Sample from rtnorm (l=",l,", u=",u,", m=",m,", sd=",s,")"),cex=cex,cex.main=cex.main)