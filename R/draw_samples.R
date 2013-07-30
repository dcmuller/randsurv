#' @title Sample survival times from a Weibull distribution.
#' 
#' @description \code{weibull_simp} simulates survial times from a Weibull 
#' distribution given a baseline hazard, a shape parameter. Optionally,
#' covariates and regression coefficients (log hazard ratios) can be supplied.
#' This function is simple in the sense that it assumes an absense of
#' competing risks. See \code{\link{weibull_compet}} if you wish to sample
#' survival times for competing events.
#' 
#' Note well that the rate and shape parameters do not correspond to those
#' in the base function \code{\link{rweibull}}. The parametrisation implemented
#' here samples from a Weibull distribution with hazard function 
#' \eqn{h(t) = p \lambda^{p-1}}{h(t) = p * lambda^(p-1)}, where 
#' \eqn{\lambda = lambda_0\exp(x\beta)}{lambda = lambda0 * exp(xb)}.
#' 
#' @param lambda0 the rate parameter for the baseline hazard
#' @param p the shape parameter for the baseline hazard
#' @param cens_time censor time beyond which failures are not observed
#' @param X optional matrix of regressors
#' @param beta optional vector or regression coefficients (log hazard ratios)
#' 
#' @export
#' 
#' @examples
#' require(survival)
#' X <- as.matrix(rbinom(100, size = 1, prob = .5))
#' beta <- c(2)
#' fupdata <- weibull_simp(lambda0 = 10e-7, p = 3, cens_time=50, X=X, beta=beta)
#' plot(survfit(Surv(fupdata$obs_time, fupdata$status) ~ X))
#' 
#' @author David C Muller
weibull_simp <- function(lambda0, p, cens_time, X = NULL , beta = NULL) {
  if (is.null(X) && !is.null(beta)) {
    stop("You have specified regression coefficients 'beta', but no 
         matrix of regressors 'X'")
    stop(paste0("You have specified regression coefficients 'beta, but", 
                "\n",
                "no matrix of regressors 'X'"
        )
    )
  }
  else if (is.null(beta) && !is.null(X)) {
    stop(paste0("You have specified a matrix of regressors 'X', but",  
                "\n",
                "no regression coefficients 'beta'"
        )
    )
  }
  else if (is.null(X) && is.null(beta)) { 
    lambda <- lambda0
  }
  else {
    lambda <- lambda0 * exp(X %*% beta)
  }
  ftime <- vector(mode="numeric", length=nrow(X))
  for (i in 1:length(ftime)) {
     ftime[i] <- inv_weibull(lambda[i], p, rand_unif=runif(1))
  }
  time <- pmin(ftime, cens_time)
  status <- as.numeric(time==ftime)
  fup_frame <- data.frame(f_time <- ftime, obs_time = time, status = status)
  return(fup_frame)
}

