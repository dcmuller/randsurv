#include <Rcpp.h>
using namespace Rcpp;
// dcmuller 20131021
// defines survival functions and functors (function objects) to
//    sample time to failure. Functors are necessary so that we can
//    solve for time using a generic function implementing 
//    the Newton-Raphson root finder.

// numerical derivative (template function)
template<typename TYPE>
double num_deriv( TYPE f,
                  double x,
                  double dx = 10e-7)
{
  double d_dx = (f(x + dx) - f(x)) / dx;
  return d_dx;
}

// Newton-Raphson algorithm (template function)
template<typename TYPE>
double newt_raph( TYPE f, 
                  double init = 1.0, 
                  double tol  = 1e-10, 
                  int maxiter = 500)
{
  double xnew, x, dfdx;
  int i = 0;
  x = init;
  while (i < maxiter) {
    dfdx = num_deriv(f, x);
    xnew = x - (f(x) / dfdx);
    if (abs(x - xnew) < tol ) break;
    i++;
    x = xnew;
  }
  if (i >= maxiter) {
    stop("exceeded maximum number of iterations");
  }
  else {
    return xnew;
  }
}

// Weibull survival function
// [[Rcpp::export]]
double surv_weibull(double t, double lambda, double p, double rand_unif) {
  return (lambda * pow(t,p)) + log(rand_unif);
}

class surv_weibull_call {
  private:
    double lambda, p, rand_unif;
    
  public:
    surv_weibull_call(double _lambda, double _p, double _rand_unif) :
      lambda(_lambda), p(_p), rand_unif(_rand_unif) {};
    
    double operator()(double t) const {
      return surv_weibull(t, lambda, p, rand_unif);
    }
};

// [[Rcpp::export]] 
double nrsolve_surv_weibull( double lambda, double p, double rand_unif, 
                      double init_t = 1, double tol = 1e-10, int maxiter = 500) 
{
  // instantiate surv_weibull functor with supplied arguments
  surv_weibull_call f(lambda, p, rand_unif);
  
  // solve for t
  double t = newt_raph(f, init_t, tol, maxiter);
  return t;
}
    
    
// derivative of Weibull survival function
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
class surv_weibull_compet_call {
  private:
    double lambda_1, p_1, lambda_2, p_2, rand_unif;
    
  public:
    surv_weibull_compet_call(
      double _lambda_1, 
      double _p_1,
      double _lambda_2,
      double _p_2,
      double _rand_unif) :
        lambda_1(_lambda_1), 
        p_1(_p_1), 
        lambda_2(_lambda_2), 
        p_2(_p_2), 
        rand_unif(_rand_unif) {};
    
    double operator()(double t) const {
      return surv_weibull_compet(t, lambda_1, p_1, lambda_2, p_2, rand_unif);
    }
};

// [[Rcpp::export]] 
double nrsolve_surv_weibull_compet( double lambda_1, double p_1, 
                                    double lambda_2, double p_2,
                                    double rand_unif, double init_t = 1, 
                                    double tol = 1e-10, int maxiter = 500) 
{
  // instantiate surv_weibull functor with supplied arguments
  surv_weibull_compet_call f(lambda_1, p_1, lambda_2, p_2, rand_unif);
  
  // solve for t
  double t = newt_raph(f, init_t, tol, maxiter);
  return t;
}

// derivative of Weibull survival function with competing risks
// [[Rcpp::export]]
double ddt_surv_weibull_compet( double t, 
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
