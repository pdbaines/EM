\name{MCEM}
\alias{MCEM}
\title{Function for executing an instance of the Monte Carlo EM algorithm}
\usage{
  MCEM(y.obs, theta.0, fixed, update, max.iter,
    monitor = TRUE, ss.MC.reps = 10, mu.start = 0,
    sigma.start = 1e+05, iter.to.wait = 5, cr = 0.95,
    n.monitor.samples = 10000, smooth = FALSE, penalty = 2,
    B = 100, iter.to.wait.sp = 50, spline.se.every = 100,
    logLike = NULL, keep.paths = TRUE,
    append.paths = FALSE, tol = 1e-10,
    tol.type = "relative", print.every = Inf,
    compute.DM = FALSE, verbose = FALSE)
}
\arguments{
  \item{y.obs}{observed data (any format)}

  \item{theta.0}{initial value for parameter}

  \item{fixed}{any additional fixed quantities (mainly for
  computational efficiency)}

  \item{update}{function that outputs theta^{(t+1)} given
  theta^{(t)}. Must take take arguments \code{theta.t},
  \code{y.obs}, \code{fixed} and \code{verbose} (in that
  order) i.e., it will be called as: \code{theta.t1 <-
  update(theta.t,y.obs,fixed,verbose)} When using
  \code{MCEM} it is expected that \code{update} will be a
  stochastic function, and thus return (possibly) different
  values even if called with the same input arguments.}

  \item{max.iter}{maximum number of iterations to run the
  algorithm}

  \item{monitor}{whether to perform advanced convergence
  monitoring using the methods in Baines, Xu and Wang
  (2013). If \code{FALSE} then basic (but potentially
  unreliable) convergence monitoring using standard
  relative or absolute tolerance is used. If \code{TRUE}
  then the estimate is taken as the fixed point of the
  update mapping, which is approximated via either a linear
  or spline fit to the EM path.}

  \item{ss.MC.reps}{the number of calls to \code{update}
  used to estimate the Monte-Carlo variability in the
  update function.}

  \item{mu.start}{starting value for mu}

  \item{sigma.start}{starting value for sigma}

  \item{iter.to.wait}{specifies the number of iterations to
  wait before checking for convergence begins using a
  linear fit}

  \item{cr}{the confidence level used in convergence
  checking}

  \item{n.monitor.samples}{the number of monitoring samples
  to use}

  \item{smooth}{whether to use a smooth function to
  approximate the update mapping. If \code{TRUE} then a
  spline fit is used, otherwise a simple linear fit is
  used.}

  \item{penalty}{the penalty term used if a smooth function
  is used to approximate the update mapping}

  \item{B}{the bootstrap sample size used for approximating
  the update mapping}

  \item{iter.to.wait.sp}{the number of iterations to wait
  before attempting to compute a spline fit to the update
  mapping}

  \item{spline.se.every}{the number of iterations at which
  the SE is recomputed from the spline fit to the update
  mapping. Since this can be slow, it is not recommended to
  compute this at each iteration.}

  \item{logLike}{(optional) function to compute the
  log-likelihood (or log-posterior) at each iteration. Must
  take only two arguments: \code{theta.t} and \code{y.obs}.
  Will be called as e.g., \code{logLike(theta.t,y.obs)}}

  \item{keep.paths}{whether to store the value of the
  parameter at each iteration, or just return the final
  converged value.}

  \item{append.paths}{only used if \code{keep.paths=TRUE}.
  if \code{TRUE} then parameter values are appended after
  each iteration. This is fast initially, but if lots of
  iterations are required it will slow down. If
  \code{FALSE}, a full block of size \code{max.iter} is
  allocated to store the parameter values, which can be
  unnecessary if the algorithm converges quickly.}

  \item{tol}{the tolerance used when checking for
  convergence at each iteration}

  \item{tol.type}{the type of stopping rule used. Options
  are \code{relative} and \code{absolute}.}

  \item{print.every}{specifies the interval at which status
  updates are printed to the screen e.g.,
  \code{print.every=2} will print to screen after every two
  iterations are completed.}

  \item{compute.DM}{now defunct}

  \item{verbose}{integer controling the level of verbosity
  of the function}
}
\value{
  Lots of stuff
}
\description{
  The function \code{\link{MCEM}} takes an \code{update}
  function, observed data \code{y.obs}, starting values
  \code{theta.0}, any fixed quantities \code{fixed} (see
  below) and other arguments and runs the algorithm until a
  stopping rule is satisfied or a maximum number of
  iterations is reached.
}
\details{
  Lots of stuff
}
\examples{
"MCEM.example.update" <- function(theta,y.obs,fixed,verbose)
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

t3 <- MCEM(y.obs=y.obs, fixed=fixed, theta.0=theta.0, update=EM.example.update,
           max.iter=max.iter, monitor=TRUE, ss.MC.reps=ss.MC.reps,
           mu.start=mu.start, sigma.start=sigma.start,
           smooth=smooth,B=B,penalty=penalty,spline.se.every=spline.se.every,
           iter.to.wait=iter.to.wait, cr=cr,n.monitor.samples=n.monitor.samples,
           print.every=print.every, tol=zero.tol, verbose=verbose)
print(t3)
}
\seealso{
  \code{\link{EM}}
}

