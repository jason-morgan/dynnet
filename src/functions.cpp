#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double llik_logit(NumericVector y, NumericVector lp)
{
  int n = y.size();
  double llik = 0;

  for (int i = 0; i < n; ++i) {
    double p = R::plogis(lp[i], 0.0, 1.0, 1.0, false);
    llik += R::dbinom(y[i], 1.0, p, true);
  }

  return(llik);
}

// [[Rcpp::export]]
double llik_poisson(NumericVector y, NumericVector lp)
{
  int n = y.size();
  double llik = 0;

  for (int i = 0; i < n; ++i) {
    double lambda = exp(lp[i]);
    llik += R::dpois(y[i], lambda, true);
  }

  return(llik);
}

// [[Rcpp::export]]
NumericVector norm_euclidean(NumericMatrix X)
{
  int nrow = X.nrow();
  int ncol = X.ncol();

  NumericVector norm(nrow);

  for (int i = 0; i < nrow; ++i) {
    double sumsqr = 0;
    for (int j = 0; j < ncol; ++j) {
      sumsqr += X(i, j) * X(i, j);
    }
    norm[i] = sqrt(sumsqr);
  }

  return norm;
}

// [[Rcpp::export]]
double distance_penalty(NumericMatrix X)
{
  int nrow = X.nrow();
  int ncol = X.ncol();

  double penalty = 0;

  for (int i = 0; i < nrow; ++i) {
    double sumsqr = 0;
    for (int j = 0; j < ncol; ++j) {
      sumsqr += X(i, j) * X(i, j);
    }
    penalty += log(sqrt(sumsqr));
  }

  return penalty;
}
