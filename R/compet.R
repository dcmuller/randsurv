# Weibull survival function
surv_weibull_compet2 <- function(t, lambda_1, p_1, lambda_2, p_2, rand_unif) {
  return((lambda_1 * t^p_1) + (lambda_2 * t^p_2) + log(rand_unif))
}

# derivative of Weibull survival function
ddt_surv_weibull_compet2 <- function(t, lambda_1, p_1, lambda_2, p_2, rand_unif) {
  return((p_1 * lambda_1 * t^(p_1-1)) + (p_2 * lambda_2 * t^(p_2-1)))
}

# which 'event' is observed is a function of the magnitudes 
# of the competing hazards
observed_event <- function(t, lambda_1, p_1, lambda_2, p_2) {
# all operations are elementwise  
  haz1 <- p_1 * lambda_1 * t^(p_1-1) 
  haz2 <- p_2 * lambda_2 * t^(p_2-1)
  pr1 <- haz1/(haz1 + haz2)
  event <- rbinom(length(t), 1, pr1)
  event[event==0] <- 2
  return(event)
}
