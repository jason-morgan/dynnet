#include <Rcpp.h>

// #include <RcppArmadillo.h>
// #include <RcppArmadilloExtensions/sample.h>

// // [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export(.C_insert_ref)]]
NumericMatrix C_insert_ref(IntegerVector ref_idx,
			 NumericMatrix ref_pos,
			 IntegerVector est_idx,
			 NumericMatrix est_pos)
{
  NumericMatrix result(ref_idx.size() + est_idx.size(),
		       ref_pos.ncol());

  // R indices start at 1
  ref_idx = ref_idx - 1;
  est_idx = est_idx - 1;

  for (int i = 0; i < ref_idx.size(); ++i) {
    result(ref_idx[i],_) = ref_pos(i,_);
  }

  for (int i = 0; i < est_idx.size(); ++i) {
    result(est_idx[i],_) = est_pos(i,_);
  }

  return(result);
}

double C_euclidean(NumericVector x, NumericVector y)
{
  double d = 0.0;

  for (int i = 0; i < x.size(); ++i) {
    d += pow(x[i] - y[i], 2.0);
  }

  return(sqrt(d));
}

// [[Rcpp::export(.C_dist_euclidean)]]
NumericVector C_dist_euclidean(NumericMatrix X)
{
  int nrow = X.nrow();

  NumericVector d(nrow * (nrow-1) / 2);
  int idx = 0;

  for (int j=0; j < (nrow-1); ++j) {
    for (int i=(j+1); i < nrow; ++i) {
      d[idx] = C_euclidean(X(j,_), X(i,_));
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
