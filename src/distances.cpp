#include "dynnet.h"

using namespace Rcpp;

inline double euclidean(NumericVector x, NumericVector y)
{
  double d = 0.0;

  for (int i = 0; i < x.size(); ++i) {
    d += pow(x[i] - y[i], 2.0);
  }

  return(sqrt(d));
}

// [[Rcpp::export(.C_dist_euclidean)]]
NumericVector dist_euclidean(NumericMatrix X)
{
  int nrow = X.nrow();

  NumericVector d(nrow * (nrow-1) / 2);
  int idx = 0;

  for (int j=0; j < (nrow-1); ++j) {
    for (int i=(j+1); i < nrow; ++i) {
      d[idx] = euclidean(X(j,_), X(i,_));
      idx++;
    }
  }

  return(d);
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
