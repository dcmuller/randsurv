\name{newt_raph}
\alias{newt_raph}
\title{Univariate root finding via the Newton-Raphson algorithm}
\usage{
  newt_raph(fun, dfun = NULL, init = 1, maxiter = 500,
    tol = 1e-10, ...)
}
\arguments{
  \item{fun}{the name of the function to be solved. The
  first argument of \code{fun} is the variable to be solved
  for.}

  \item{dfun}{the name of the function that is the
  analytical derivative of \code{fun}.}

  \item{init}{initial value}

  \item{maxiter}{maximum number of Newton-Raphson
  iterations}

  \item{tol}{numerical tolerance for the solution}

  \item{...}{other arguments passed to \code{fun}}
}
\description{
  This function applies the Newton-Raphson algorithm to
  find the root of a given function. Numerical derivatives
  are used by default, but optionally analytical
  derivatives can be supplied.
}
\examples{
f <- function(t,lambda1,p1,lambda2,p2,rand) {
 lambda1*(t^p1)  + lambda2*(t^p2) + log(rand)
}
fprime <- function(t,lambda1,p1,lambda2,p2, rand) {
  p1*lambda1*(t^(p1-1))  + p2*lambda2*(t^(p2-1))
}
newt_raph(f, fprime, p1=1, p2=1, lambda1=1, lambda2=1, rand=.2)
newt_raph(f, p1=1, p2=1, lambda1=1, lambda2=1, rand=.2)
}

