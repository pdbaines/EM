\name{EM}
\alias{EM}
\title{Function for executing an instance of the EM algorithm}
\usage{
  EM(y.obs, theta.0, fixed, update, max.iter,
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
  update(theta.t,y.obs,fixed,verbose)}}

  \item{max.iter}{maximum number of iterations to run the
  algorithm}

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
  The function \code{\link{EM}} takes an \code{update}
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
# From (ref)
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
max.iter <- 1000
verbose <- FALSE
t1 <- EM(y.obs=y.obs, fixed=fixed, theta.0=theta.0, update=EM.example.update, max.iter=max.iter, verbose=verbose)
t1
}
\seealso{
  \code{\link{MCEM}}
}

