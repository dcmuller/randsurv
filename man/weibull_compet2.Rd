\name{weibull_compet2}
\alias{weibull_compet2}
\title{Sample event times from a Weibull distribution with competing risks}
\usage{
  weibull_compet2(lambda0, p, cens_time, X_1 = NULL,
    beta_1 = NULL, X_2 = NULL, beta_2 = NULL)
}
\arguments{
  \item{lambda0}{vector of length 2 containing the rate
  parameters for the baseline sub-hazards}

  \item{p}{vector of length 2 containing the shape
  parameters for the baseline sub-hazards}

  \item{cens_time}{censor time beyond which events are not
  observed}

  \item{X_1}{optional matrix of regressors for the first
  competing event}

  \item{beta_1}{optional vector of regression coefficients
  (log hazard ratios) for the first competing event}

  \item{X_2}{optional matrix of regressors for the second
  competing event}

  \item{beta_2}{optional vector of regression coefficients
  (log hazard ratios) for the second competing event}
}
\description{
  \code{weibull_compet2} simulates event times for two
  competing events given Weibull sub-hazards. Optionally,
  covariates and regression coefficients (log hazard
  ratios) can be supplied for either or both of the
  competing events. See \code{\link{weibull_simp}} if you
  wish to sample survival times for one event in the
  absence of competing risks.

  Note well that the rate and shape parameters of the
  Weibull sub-hazards do not correspond to those in the
  base function \code{\link{rweibull}}. The parametrisation
  implemented here samples from a Weibull distribution with
  sub-hazard function \eqn{h(t) = p \lambda^{p-1}}{h(t) = p
  * lambda^(p-1)}, where \eqn{\lambda =
  lambda_0\exp(x\beta)}{lambda = lambda0 * exp(xb)}.
}
\examples{
require(survival)
X1 <- as.matrix(rbinom(1000, size = 1, prob = .5))
X2 <- as.matrix(rbinom(1000, size = 1, prob = .1))
beta1 <- c(log(2))
beta2 <- c(log(3))
fupdata <- weibull_compet2(lambda0 = c(10e-10, 10e-9), p = c(4.3, 4.2),
   cens_time=80, X_1 = X1, beta_1 = beta1, X_2 = X2, beta_2 = beta2)
plot(survfit(Surv(fupdata$obs_time, fupdata$status==1) ~ X1))
}
\author{
  David C Muller
}

