\name{weibull_simp}
\alias{weibull_simp}
\title{Sample survival times from a Weibull distribution.}
\usage{
  weibull_simp(lambda0, p, cens_time, X = NULL,
    beta = NULL)
}
\arguments{
  \item{lambda0}{the rate parameter for the baseline
  hazard}

  \item{p}{the shape parameter for the baseline hazard}

  \item{cens_time}{censor time beyond which failures are
  not observed}

  \item{X}{optional matrix of regressors}

  \item{beta}{optional vector or regression coefficients
  (log hazard ratios)}
}
\description{
  \code{weibull_simp} simulates survial times from a
  Weibull distribution given a baseline hazard, a shape
  parameter. Optionally, covariates and regression
  coefficients (log hazard ratios) can be supplied. This
  function is simple in the sense that it assumes an
  absense of competing risks. See
  \code{\link{weibull_compet}} if you wish to sample
  survival times for competing events.

  Note well that the rate and shape parameters do not
  correspond to those in the base function
  \code{\link{rweibull}}. The parametrisation implemented
  here samples from a Weibull distribution with hazard
  function \eqn{h(t) = p \lambda^{p-1}}{h(t) = p *
  lambda^(p-1)}, where \eqn{\lambda =
  lambda_0\exp(x\beta)}{lambda = lambda0 * exp(xb)}.
}
\examples{
require(survival)
X <- as.matrix(rbinom(100, size = 1, prob = .5))
beta <- c(2)
fupdata <- weibull_simp(lambda0 = 10e-7, p = 3, cens_time=50, X=X, beta=beta)
plot(survfit(Surv(fupdata$obs_time, fupdata$status) ~ X))
}
\author{
  David C Muller
}

