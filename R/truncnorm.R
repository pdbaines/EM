
##
## Truncated (univariate) normal distribution functions
## (Note: R has library(msm) containing equivalent functions...)
##

#' @export
"dtnorm" <- function(x,lower=-Inf,upper=Inf,mean=0,sd=1,log=FALSE)
{
  num <- dnorm(x=x,mean=mean,sd=sd,log=log)
  den <- pnorm(q=upper,mean=mean,sd=sd)-pnorm(q=lower,mean=mean,sd=sd)
  if (log){
    ret <- num - log(den)
  } else {
    ret <- num/den
  }
  # Take care of out of range cases:
  ret[x<lower | x>upper] <- ifelse(log,-Inf,0)
  return(ret)
}

#' @export
"rtnorm" <- function(n,lower=-Inf,upper=Inf,mean=0,sd=1)
{
  return(qnorm(runif(n=n,min=pnorm(lower),max=pnorm(upper))))
}

#' @export
"ptnorm" <- function(q,lower=-Inf,upper=Inf,mean=0,sd=1,lower.tail=TRUE,log.p=FALSE)
{
  # Note this may not be numerically stable for very small/large percentiles...
  val <- pnorm(q=q,mean=mean,sd=sd)
  low <- pnorm(q=lower,mean=mean,sd=sd)
  upp <- pnorm(q=upper,mean=mean,sd=sd)
  ret <- (val-low)/(upp-low)
  ret[q<lower] <- 0.0
  ret[q>upper] <- 1.0
  if (log.p)
    ret <- log(ret)
  return(ret)
}

#' @export
"qtnorm" <- function(p,lower=-Inf,upper=Inf,mean=0,sd=1,lower.tail=TRUE,log.p=FALSE)
{
  upp <- pnorm(q=upper,mean=mean,sd)
  low <- pnorm(q=lower,mean=mean,sd)
  if (log.p){
    val <- low + exp(p)*(upp-low)
  } else {
    val <- low + p*(upp-low)
  }
  ret <- qnorm(p=val,mean=mean,sd=sd)
  return(ret)
}

