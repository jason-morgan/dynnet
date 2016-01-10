#ifndef __DYNNET_H__
#define __DYNNET_H__

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

// likelihoods
double C_llik_logit(NumericVector y, NumericVector lp);
double C_llik_poisson(NumericVector y, NumericVector lp);

// distributions
arma::vec C_dmvnorm(arma::mat x,
		    arma::rowvec mean,
		    arma::mat sigma,
		    bool logd = false);

#endif
