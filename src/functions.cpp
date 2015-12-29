#include <Rcpp.h>
using namespace Rcpp;


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
