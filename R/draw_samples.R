#' @title Sample survival times from a Weibull distribution.
#' 
#' @description \code{weibull_simp} simulates survial times from a Weibull 
#' distribution given a baseline hazard, a shape parameter. Optionally,
#' covariates and regression coefficients (log hazard ratios) can be supplied.
#' This function is simple in the sense that it assumes an absense of
#' competing risks. See \code{\link{weibull_compet}} if you wish to sample
#' survival times for two competing events.
#' 
#' Note well that the rate and shape parameters do not correspond to those
#' in the base function \code{\link{rweibull}}. The parametrisation implemented
#' here samples from a Weibull distribution with hazard function 
#' \eqn{h(t) = p \lambda^{p-1}}{h(t) = p * \lambda^(p-1)}, where 
#' \eqn{\lambda = \lambda_0\exp(x\beta)}{lambda = lambda0 * exp(xb)}.
#' 
#' @param lambda0 the rate parameter for the baseline hazard
#' @param p the shape parameter for the baseline hazard
#' @param cens_time censor time beyond which failures are not observed
#' @param X optional matrix of regressors
#' @param beta optional vector of regression coefficients (log hazard ratios)
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
  fup_frame <- data.frame(f_time = ftime, obs_time = time, status = status)
  return(fup_frame)
}


#' @title Sample event times from a Weibull distribution with competing risks
#' 
#' @description \code{weibull_compet} simulates event times for two competing
#' events given Weibull sub-hazards. Optionally, covariates and regression 
#' coefficients (log hazard ratios) can be supplied for either or both of the
#' competing events. See \code{\link{weibull_simp}} if you wish to sample
#' survival times for one event in the absence of competing risks.
#' 
#' Note well that the rate and shape parameters of the Weibull sub-hazards do 
#' not correspond to those in the base function \code{\link{rweibull}}. The 
#' parametrisation implemented here samples from a Weibull distribution with 
#' sub-hazard function \eqn{h(t) = p \lambda^{p-1}}{h(t) = p * \lambda^(p-1)}, 
#' where \eqn{\lambda = \lambda_0\exp(x\beta)}{lambda = lambda0 * exp(xb)}.
#' 
#' @param lambda0 vector of length 2 containing the rate parameters for the 
#' baseline sub-hazards
#' @param p vector of length 2 containing the shape parameters for the baseline 
#' sub-hazards
#' @param cens_time censor time beyond which events are not observed
#' @param X_1 optional matrix of regressors for the first competing event
#' @param beta_1 optional vector of regression coefficients (log hazard ratios)
#' for the first competing event
#' @param X_2 optional matrix of regressors for the second competing event
#' @param beta_2 optional vector of regression coefficients (log hazard ratios)
#' for the second competing event
#' 
#' @export
#' 
#' @examples
#' require(survival)
#' X1 <- as.matrix(rbinom(1000, size = 1, prob = .5))
#' X2 <- as.matrix(rbinom(1000, size = 1, prob = .1))
#' beta1 <- c(log(2))
#' beta2 <- c(log(3))
#' fupdata <- weibull_compet(lambda0 = c(10e-10, 10e-9), p = c(4.3, 4.2), 
#'    cens_time=80, X_1 = X1, beta_1 = beta1, X_2 = X2, beta_2 = beta2)
#' plot(survfit(Surv(fupdata$obs_time, fupdata$status==1) ~ X1))
#' 
#' @author David C Muller
weibull_compet <- function(lambda0, p, cens_time, 
                            X_1 = NULL, beta_1 = NULL,
                            X_2 = NULL, beta_2 = NULL
                            ) {
  if (is.null(X_1) && !is.null(beta_1)) {
    stop(paste0("You have specified regression coefficients 'beta_1', but", 
                "\n",
                "no matrix of regressors 'X_1'"
        )
    )
  }
  else if (is.null(beta_1) && !is.null(X_1)) {
    stop(paste0("You have specified a matrix of regressors 'X_1', but",  
                "\n",
                "no regression coefficients 'beta_1'"
        )
    )
  }
  else if (is.null(X_1) && is.null(beta_1)) { 
    lambda_1 <- lambda0[1]
  }
  else {
    lambda_1 <- lambda0[1] * exp(X_1 %*% beta_1)
  }
  if (is.null(X_2) && !is.null(beta_2)) {
    stop(paste0("You have specified regression coefficients 'beta_2', but", 
                "\n",
                "no matrix of regressors 'X_1'"
        )
    )
  }
  else if (is.null(beta_2) && !is.null(X_2)) {
    stop(paste0("You have specified a matrix of regressors 'X_2', but",  
                "\n",
                "no regression coefficients 'beta_2'"
        )
    )
  }
  else if (is.null(X_2) && is.null(beta_2)) { 
    lambda_2 <- lambda0[2]
  }
  else {
    lambda_2 <- lambda0[2] * exp(X_2 %*% beta_2)
  }
  ftime <- vector(mode="numeric", length=nrow(X_1))
  event <- vector(mode="numeric", length=nrow(X_1))
  for (i in 1:length(ftime)) {
    ftime[i] <- newt_raph(randsurv:::surv_weibull_compet, # S(t) 
                          randsurv:::ddt_surv_weibull_compet, # dS(t)/dt
                          lambda_1 = lambda_1[i],
                          p_1 = p[1],
                          lambda_2 = lambda_2[i],
                          p_2 = p[2], 
                          rand_unif = runif(1)
    )
    event[i] <- observed_event(t = ftime[i], 
                               lambda_1 = lambda_1[i], 
                               p_1 = p[1],
                               lambda_2 = lambda_2[i], 
                               p_2 = p[2]
    ) 
  }
  time <- pmin(ftime, cens_time)
  status <- as.numeric(time==ftime)*event
  fup_frame <- data.frame(f_time = ftime, obs_time = time, status = status)
  return(fup_frame)
}
