# invert CDF to draw survival time
# Weibull survival function
surv_weibull <- function(t, lambda, p, rand_unif) {
  return((lambda * t^p) + log(rand_unif))
}
# derivative of Weibull survival function
ddt_surv_weibull <- function(t,lambda, p, rand_unif) {
  return(p * lambda * t^(p-1))
}
# use explicit inversion of Weibull survival function to solve for time
inv_weibull <- function(lambda, p, rand_unif) {
  time <- (-log(rand_unif)/lambda)^(1/p)
  return(time)
}
