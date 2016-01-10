#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "dynnet.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export(.C_propose_Z)]]
NumericMatrix C_propose_Z(NumericMatrix Z)
{
  int n = 1;
  int R = Z.nrow();
  int C = Z.ncol();
  NumericMatrix result(clone(Z));
  IntegerVector row_idx = Range(0, R-1);
  IntegerVector idx = Rcpp::RcppArmadillo::sample(row_idx, n, FALSE,
						  NumericVector::create());
  for (int i = 0; i < n; ++i) {
    result(idx[i],_) = Z(idx[i],_) + rnorm(C, 0.0, 1.0/R);
  }

  return(result);
}

// [[Rcpp::export(.C_propose_beta0)]]
double C_propose_beta0(double beta0)
{
  return(std::abs(beta0 + R::rnorm(0.0, 0.8)));
}

// [[Rcpp::export(.C_log_prior_beta0)]]
double C_log_prior_beta0(double beta0)
{
  // Matches the prior used by HRH (2002)
  return(R::dgamma(beta0, 1.0, 1.0, 1));
}

// [[Rcpp::export(.C_log_prior_Z)]]
double C_log_prior_Z(NumericMatrix Z)
{
  int n = Z.ncol();
  NumericVector mu(n);
  NumericMatrix sigma(n);

  mu.fill(0.0);
  sigma.fill_diag(2.0);

  return(sum(C_dmvnorm(as<arma::mat>(Z), as<arma::rowvec>(mu),
		       as<arma::mat>(sigma), true)));
}

// [[Rcpp::export(.C_log_posterior_logit)]]
double C_log_posterior_logit(NumericVector y, NumericVector lp,
			     double beta0, NumericMatrix Z)
{
  double post = C_llik_logit(y, lp)
    + C_log_prior_beta0(beta0) + C_log_prior_Z(Z);

  return(post);
}
