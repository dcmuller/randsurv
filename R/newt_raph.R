#' @title Numerical derivative of a function
#' 
#' @description Returns df/dx evaluated at \code{x} for a given function 
#' \code{f(...)}.
#'
#' @param fun the function to numerically differentiate. The derivative is
#' calculated with respect to the first argument of \code{fun}.
#' @param x the value at which to evaluate the derivative
#' @param dx the width of the differential, defaults to \code{10e-9}.
#' @param ... other arguments passed to \code{fun}
num_deriv <- function(fun, x, dx=10e-9, ...) {
  d_dx <- (fun(x + dx, ...) - fun(x, ...)) / dx
  return(d_dx)
}


#' @title Univariate root finding via the Newton-Raphson algorithm
#' 
#' @description This function applies the Newton-Raphson algorithm to find the 
#' root of a given function. Numerical derivatives are used by default, but
#' optionally analytical derivatives can be supplied.
#' 
#' @param fun the name of the function to be solved. The first argument of 
#' \code{fun} is the variable to be solved for.
#' @param dfun the name of the function that is the analytical derivative of
#' \code{fun}.
#' @param init initial value
#' @param maxiter maximum number of Newton-Raphson iterations
#' @param tol numerical tolerance for the solution
#' @param ... other arguments passed to \code{fun}
#' @export
#' @examples
#' f <- function(t,lambda1,p1,lambda2,p2,rand) {
#'  lambda1*(t^p1)  + lambda2*(t^p2) + log(rand)
#' }
#' fprime <- function(t,lambda1,p1,lambda2,p2, rand) {
#'   p1*lambda1*(t^(p1-1))  + p2*lambda2*(t^(p2-1))
#' }
#' newt_raph(f, fprime, p1=1, p2=1, lambda1=1, lambda2=1, rand=.2)
#' newt_raph(f, p1=1, p2=1, lambda1=1, lambda2=1, rand=.2)

newt_raph <- function(fun,
                      dfun=NULL,
                      init=1,
                      maxiter=500,
                      tol=1e-10, 
                      ...
                      ) {
  if (is.null(dfun)) {
    # no analytical derivative specified, use numerical differentiation
    dfun <- function(...) num_deriv(fun=fun, ...)
  }
  i <- 1
  x <- init
  while (i <= maxiter) {
    xnew <- x - fun(x,...)/dfun(x, ...)
    if (abs(x - xnew) < tol) break
    i <- i+1
    x <- xnew
  }
  if (i >= maxiter) {
    stop("exceeded maximum number of iterations")
  }
  else {
    return(xnew)
  }
}
