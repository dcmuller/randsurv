#include <Rcpp.h>
using namespace Rcpp;

// Weibull survival function
// [[Rcpp::export]]
double surv_weibull(double t, double lambda, double p, double rand_unif) {
  return (lambda * pow(t,p)) + log(rand_unif);
}

// derivative of Weibull survival function
// [[Rcpp::export]]
double ddt_surv_weibull (double t, double lambda, double p, double rand_unif) {
  return p * lambda * pow(t, p-1);
}

// use explicit inversion of Weibull survival function to solve for time
// [[Rcpp::export]]
double inv_weibull(double lambda, double p, double rand_unif) {
  double time;
  time = pow((-log(rand_unif)/lambda), (1/p));
  return time;
}

// Weibull survival function with competing risks
// [[Rcpp::export]]
double surv_weibull_compet(double t, 
                           double lambda_1, 
                           double p_1, 
                           double lambda_2, 
                           double p_2, 
                           double rand_unif) {
  return (lambda_1 * pow(t, p_1)) + (lambda_2 * pow(t, p_2)) + log(rand_unif);
}

// derivative of Weibull survival function with competing risks
// [[Rcpp::export]]
double ddt_surv_weibull_compet(double t, 
                                double lambda_1, 
                                double p_1, 
                                double lambda_2, 
                                double p_2, 
                                double rand_unif) {
  return (p_1 * lambda_1 * pow(t, (p_1-1))) + (p_2 * lambda_2 * pow(t, (p_2-1)));
}

// which 'event' is observed is a function of the magnitudes 
// of the competing hazards
// [[Rcpp::export]]
int observed_event( double t, 
                    double lambda_1, 
                    double p_1, 
                    double lambda_2, 
                    double p_2) {
  RNGScope scope;
  int fail, event;
  double haz1, haz2, pr1;
  haz1 = p_1 * lambda_1 * pow(t, (p_1-1));
  haz2 = p_2 * lambda_2 * pow(t, (p_2-1));
  pr1 = haz1/(haz1 + haz2);
  fail = int(R::rbinom(1, pr1));
  event = fail==1 ? 1 : 2;
  return event;
}
