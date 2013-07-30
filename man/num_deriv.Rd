\name{num_deriv}
\alias{num_deriv}
\title{Numerical derivative of a function}
\usage{
  num_deriv(fun, x, dx = 1e-08, ...)
}
\arguments{
  \item{fun}{the function to numerically differentiate. The
  derivative is calculated with respect to the first
  argument of \code{fun}.}

  \item{x}{the value at which to evaluate the derivative}

  \item{dx}{the width of the differential, defaults to
  \code{10e-9}.}

  \item{...}{other arguments passed to \code{fun}}
}
\description{
  Returns df/dx evaluated at \code{x} for a given function
  \code{f(...)}.
}

