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
	      String family,
	      String dist_metric)
{
  LSMModel Model = {y, X, Z_idx, k, d, burnin, samplesize, interval};
  if (family == "bernoulli") {
    Model.lsm_posterior_fn = log_posterior_logit;
    Model.lsm_lpr_fn = llik_logit;
  }

  // Set the distance function
  if (dist_metric == "euclidean") {
    Model.lsm_dist_fn = dist_euclidean;
  }
  else if (dist_metric == "euclidean2") {
    Model.lsm_dist_fn = dist_euclidean2;
  }

  // initiate starting values and State
  NumericVector Xb = wrap(as<arma::mat>(X) * as<arma::vec>(beta));
  NumericVector dist = Model.lsm_dist_fn(Z);
  NumericVector lp = (Xb - dist);
  double lpr_graph = Model.lsm_lpr_fn(&Model, lp);

  LSMState State = {
    beta,
    beta,			// old_beta
    Z,
    Z,				// old_Z
    Xb,
    Xb,				// old_Xb
    dist,
    dist,			// old_dist
    lpr_graph,
    0.0,			// posterior
    0.0,			// beta_log_prior
    0,				// beta_accept
    0.0,			// Z_log_prior
    0,				// Z_accept
    1.0				// Z_proposal_sd
  };

  State.posterior = Model.lsm_posterior_fn(&Model, &State);
  State.beta_log_prior = log_prior_beta(&Model, &State);
  State.Z_log_prior = log_prior_Z(&Model, &State);

  // matrix for samples. the extra column is for storing the log probability of
  // the graph
  NumericMatrix samples(Model.samplesize,
			1 + Model.k + ((State.Z).nrow() * (State.Z).ncol()));

  // burnin
  int Z_accept_last = 0;	// records the last Z acceptance
  int niter = 2000;
  Rcpp::Rcout << std::endl << "=========="
	      << std::endl << "  Burnin"
	      << std::endl << "==========" << std::endl;
  for (int i=0; i < Model.burnin; ++i) {
    // The first check here is just to prevent a floating point exception when
    // the burnin is very small. This is a hack and should be fixed in a more
    // elegant way.
    if ((Model.burnin > 1000) && (i % (Model.burnin / 10) == 0)) {
      Rcpp::checkUserInterrupt();
      msg_mcmc_iter(i, Model.burnin, State.beta_accept, State.Z_accept);
    }

    // Adjust Z proposal SD
    if ((Model.burnin > niter) && (i % niter == 0)) {
      adapt_Z_proposal_sd(&State, niter, Z_accept_last);
      Z_accept_last = State.Z_accept;
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
    // As above, this first check here is just to prevent a floating point
    // exception when the samplesize is very small.
    if ((Model.samplesize > 1000) && (s % (Model.samplesize / 10) == 0)) {
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
  int n = 1 +
    (State->beta).size() +
    ((State->Z).nrow() * (State->Z).ncol());
  NumericVector x(n);

  x[0] = State->lpr_graph;

  for (int j = 0; j < (State->beta).size(); ++j) {
    x[j+1] = State->beta[j];	// + 1 for the log prob
  }

  int i = (State->beta).size() + 1;
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
  NumericVector mu(Model->k);
  NumericMatrix sigma(Model->k);
  mu.fill(0.0);
  sigma.fill_diag((1.0 / Model->k) * (1.0 / Model->k));
  NumericMatrix delta = wrap(rmvnorm(1, as<arma::vec>(mu),
				     as<arma::mat>(sigma)));

  State->old_beta = clone(State->beta);
  State->old_Xb = clone(State->Xb);
  State->beta = State->beta + delta(0, _);
  State->Xb = wrap(as<arma::mat>(Model->X) * as<arma::vec>(State->beta));

  NumericVector lp = (State->Xb - State->dist);
  double new_lpr_graph = (Model->lsm_lpr_fn)(Model, lp);
  double new_beta_log_prior = log_prior_beta(Model, State);

  double p = exp((new_lpr_graph + new_beta_log_prior) -
		 (State->lpr_graph + State->beta_log_prior));
  double r = R::runif(0,1);

  if (r < p) {
    State->lpr_graph = new_lpr_graph;
    State->beta_log_prior = new_beta_log_prior;
    State->beta_accept++;
  } else {
    State->beta = State->old_beta;
    State->Xb = State->old_Xb;
  }
}

void lsm_update_Z(LSMModel *Model, LSMState *State)
{
  int R = (State->Z).nrow();
  int C = (State->Z).ncol();

  NumericVector mu(C);
  NumericMatrix sigma(C);
  mu.fill(0.0);
  sigma.fill_diag((State->Z_proposal_sd)*(State->Z_proposal_sd));
  arma::mat delta = rmvnorm(R, as<arma::vec>(mu), as<arma::mat>(sigma));

  State->old_Z = clone(State->Z);
  State->old_dist = clone(State->dist);
  State->Z = wrap(as<arma::mat>(State->Z) + delta);
  State->dist = (Model->lsm_dist_fn)(State->Z);

  NumericVector lp = (State->Xb - State->dist);
  double new_lpr_graph = (Model->lsm_lpr_fn)(Model, lp);
  double new_Z_log_prior = log_prior_Z(Model, State);

  double p = exp((new_lpr_graph + new_Z_log_prior) -
		 (State->lpr_graph + State->Z_log_prior));
  double r = R::runif(0,1);

  if (r < p) {
    State->lpr_graph = new_lpr_graph;
    State->Z_log_prior = new_Z_log_prior;
    State->Z_accept++;
  } else {
    State->Z = State->old_Z;
    State->dist = State->old_dist;
  }
}

double log_posterior_logit(LSMModel *Model, LSMState *State)
{
  NumericVector Xb = wrap(as<arma::mat>(Model->X) * as<arma::vec>(State->beta));
  NumericVector lp = (Xb - (Model->lsm_dist_fn)(State->Z));

  double post = llik_logit(Model, lp)
    + log_prior_beta(Model, State)
    + log_prior_Z(Model, State);

  return post;
}
