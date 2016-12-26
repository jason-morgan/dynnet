/* ----------------------------------------------------------------------------
** This file is part of dynnet
**
** Copyright (C) 2016 Jason W. Morgan <jason.w.morgan@gmail.com>
**
** dynnet and is free software: you can redistribute it and/or modify it under
** the terms of the GNU General Public License as published by the Free Software
** Foundation, either version 3 of the License, or (at your option) any later
** version.
**
** This program is distributed in the hope that it will be useful, but WITHOUT
** ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
** FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
** details.
**
** You should have received a copy of the GNU General Public License along with
** this program.  If not, see <http://www.gnu.org/licenses/>.
**
** -------------------------------------------------------------------------- */


#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "dynnet.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export(.C_lsm_MH)]]
List C_lsm_MH(NumericVector y,
	      NumericMatrix X,	// model matrix
	      NumericVector Z_idx,
	      int k,
	      int d,
	      int burnin,
	      int samplesize,
	      int interval,
	      NumericVector beta,
	      NumericMatrix Z,
	      String family)
{
  LSMModel Model = {y, X, Z_idx, k, d, burnin, samplesize, interval};
  if (family == "bernoulli")
    Model.lsm_posterior_fn = log_posterior_logit;

  LSMState State = {beta, Z, 0.0, 0, 0};
  NumericMatrix samples(Model.samplesize,
			Model.k + ((State.Z).nrow() * (State.Z).ncol()));
  State.posterior = Model.lsm_posterior_fn(&Model, &State);

  // burnin
  Rcpp::Rcout << std::endl << "=========="
	      << std::endl << "  Burnin"
	      << std::endl << "==========" << std::endl;
  for (int i=0; i < Model.burnin; ++i) {
    if (i % (Model.burnin / 10) == 0) {
      Rcpp::checkUserInterrupt();
      msg_mcmc_iter(i, Model.burnin, State.beta_accept, State.Z_accept);
    }

    lsm_update_Z(&Model, &State);
    lsm_update_beta(&Model, &State);
  }

  // sampling
  Rcpp::Rcout << std::endl << "=========="
	      << std::endl << " Sampling"
	      << std::endl << "==========" << std::endl;

  State.beta_accept = 0;
  State.Z_accept = 0;
  for (int s=0; s < Model.samplesize; ++s) {
    if (s % (Model.samplesize / 10) == 0) {
      Rcpp::checkUserInterrupt();
      msg_mcmc_iter(s*Model.interval, Model.samplesize*Model.interval,
		    State.beta_accept, State.Z_accept);
    }

    for (int i=0; i < Model.interval; ++i) {
      lsm_update_Z(&Model, &State);
      lsm_update_beta(&Model, &State);
    }

    save_sample(&State, &samples, s);
  }

  return Rcpp::List::create(Rcpp::Named("beta") = State.beta,
			    Rcpp::Named("Z") = State.Z,
			    Rcpp::Named("posterior") = State.posterior,
			    Rcpp::Named("beta_accept") = State.beta_accept,
			    Rcpp::Named("Z_accept") = State.Z_accept,
			    Rcpp::Named("samples") = samples);
}

void save_sample(LSMState *State, NumericMatrix *samples, int s)
{
  int n = (State->beta).size() + ((State->Z).nrow() * (State->Z).ncol());
  NumericVector x(n);

  for (int j = 0; j < (State->beta).size(); ++j) {
    x[j] = State->beta[j];
  }

  int i = (State->beta).size();
  for (int c = 0; c < (State->Z).ncol(); ++c) {
    for (int r = 0; r < (State->Z).nrow(); ++r) {
      x[i] = State->Z(r,c);
      ++i;
    }
  }

  (*samples)(s,_) = x;
}

void lsm_update_beta(LSMModel *Model, LSMState *State)
{
  NumericVector orig_beta = clone(State->beta);

  if ((State->beta).size() == 1) {
    // Keep it positive when no other covariates are included
    (State->beta)(0) = std::abs((State->beta)(0) + R::rnorm(0.0, 1.1));
  }
  else {
    State->beta = State->beta + rnorm(Model->k, 0.0, 0.05/Model->k);
  }

  double new_posterior = (Model->lsm_posterior_fn)(Model, State);

  double r = R::runif(0,1);
  double p = exp(new_posterior - State->posterior);

  if (r < p) {
    State->posterior = new_posterior;
    State->beta_accept++;
  }
  else {
    State->beta = orig_beta;
  }
}

void lsm_update_Z(LSMModel *Model, LSMState *State)
{
  int R = (State->Z).nrow();
  int C = (State->Z).ncol();
  NumericVector idx = Model->Z_idx - 1;	// R indices start at 1

  NumericMatrix orig_Z = clone(State->Z);

  for (int i=0; i < idx.size(); ++i) {
    (State->Z)(idx[i], _) = (State->Z)(idx[i], _)
      + rnorm(C, 0.0, 2.0/idx.size());
  }

  double new_posterior = (Model->lsm_posterior_fn)(Model, State);
  double r = R::runif(0, 1);
  double p = exp(new_posterior - State->posterior);

  if (r < p) {
    State->posterior = new_posterior;
    State->Z_accept++;
  }
  else {
    State->Z = orig_Z;
  }
}

double log_posterior_logit(LSMModel *Model, LSMState *State)
{
  NumericVector Xb = wrap(as<arma::mat>(Model->X) * as<arma::vec>(State->beta));
  NumericVector lp = (Xb - dist_euclidean(State->Z));

  double post = llik_logit(Model, lp)
    + log_prior_beta(Model, State)
    + log_prior_Z(Model, State);

  return post;
}

double log_prior_beta(LSMModel *Model, LSMState *State)
{
  double p = 0.0;

  if (Model->k == 1) {
    // Matches the prior used by HRH (2002)
    p = R::dgamma(State->beta(0), 1.0, 1.0, 1);
  }
  else {
    int n = Model->k;
    NumericMatrix beta(1, Model->k);
    NumericVector mu(n);
    NumericMatrix sigma(n);

    beta(0,_) = State->beta;
    mu.fill(0.0);
    sigma.fill_diag(2.0);

    p = sum(dmvnorm(as<arma::mat>(beta),
		    as<arma::rowvec>(mu),
		    as<arma::mat>(sigma), true));
  }

  return(p);
}

double log_prior_Z(LSMModel *Model, LSMState *State)
{
  int n = Model->d;
  NumericVector mu(n);
  NumericMatrix sigma(n);

  mu.fill(0.0);
  sigma.fill_diag(2.0);

  return(sum(dmvnorm(as<arma::mat>(State->Z),
		     as<arma::rowvec>(mu),
		     as<arma::mat>(sigma), true)));
}
