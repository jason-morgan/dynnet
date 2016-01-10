#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export(C_llik_logit)]]
double C_llik_logit(NumericVector y, NumericVector lp)
{
  int n = y.size();
  double llik = 0;

  for (int i = 0; i < n; ++i) {
    double p = R::plogis(lp[i], 0.0, 1.0, 1.0, false);
    llik += R::dbinom(y[i], 1.0, p, true);
  }

  return(llik);
}

// [[Rcpp::export(C_llik_poisson)]]
double C_llik_poisson(NumericVector y, NumericVector lp)
{
  int n = y.size();
  double llik = 0;

  for (int i = 0; i < n; ++i) {
    double lambda = exp(lp[i]);
    llik += R::dpois(y[i], lambda, true);
  }

  return(llik);
}
