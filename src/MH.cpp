#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "dynnet.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

LSMState copy_lsm_state(LSMState *State)
{
  LSMState x = { State->alpha,
		 State->beta,
		 State->Z,
		 State->posterior,
		 State->alpha_accept,
		 State->beta_accept,
		 State->Z_accept };
  return x;
}

// [[Rcpp::export(.C_lsm_MH)]]
List C_lsm_MH(NumericVector y,
	      NumericMatrix X,	// model matrix
	      NumericVector Z_idx,
	      int k,
	      int burnin,
	      int samplesize,
	      int interval,
	      double alpha,
	      NumericVector beta,
	      NumericMatrix Z,
	      String family)
{
  LSMModel Model = {y, X, Z_idx, k, burnin, samplesize, interval};
  if (family == "bernoulli")
    Model.lsm_posterior_fn = log_posterior_logit;

  LSMState State = {alpha, beta, Z, 0.0, 0, 0, 0};
  NumericMatrix samples(Model.samplesize,
			(1 + ((State.Z).nrow() * (State.Z).ncol())));
  State.posterior = Model.lsm_posterior_fn(&Model, &State);

  // burnin
  for (int i=0; i < Model.burnin; ++i) {
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    lsm_update_Z(&Model, &State);
    lsm_update_alpha(&Model, &State);
  }

  // sampling
  for (int s=0; s < Model.samplesize; ++s) {
    for (int i=0; i < Model.interval; ++i) {
      if (i % 1000 == 0)
	Rcpp::checkUserInterrupt();

      lsm_update_Z(&Model, &State);
      lsm_update_alpha(&Model, &State);
    }
    save_sample(&State, &samples, s);
  }

  return Rcpp::List::create(Rcpp::Named("alpha") = State.alpha,
			    Rcpp::Named("Z") = State.Z,
			    Rcpp::Named("posterior") = State.posterior,
			    Rcpp::Named("alpha_accept") = State.alpha_accept,
			    Rcpp::Named("Z_accept") = State.Z_accept,
			    Rcpp::Named("samples") = samples);
}

void save_sample(LSMState *State, NumericMatrix *samples, int s)
{
  int n = 1 + ((State->Z).nrow() * (State->Z).ncol());
  NumericVector x(n);
  x[0] = State->alpha;

  int i = 1;
  for (int c = 0; c < (State->Z).ncol(); ++c) {
    for (int r = 0; r < (State->Z).nrow(); ++r) {
      x[i] = State->Z(r,c);
      ++i;
    }
  }

  (*samples)(s,_) = x;
}

void lsm_update_alpha(LSMModel *Model, LSMState *State)
{
  LSMState Proposal = copy_lsm_state(State);
  double orig_alpha = Proposal.alpha;
  Proposal.alpha = std::abs(Proposal.alpha + R::rnorm(0.0, 1.0));
  double new_posterior = (Model->lsm_posterior_fn)(Model, &Proposal);

  double r = R::runif(0,1);
  double p = exp(new_posterior - Proposal.posterior);

  if (r < p) {
    Proposal.posterior = new_posterior;
    Proposal.alpha_accept++;
  }
  else {
    Proposal.alpha = orig_alpha;
  }

  State->alpha = Proposal.alpha;
  State->posterior = Proposal.posterior;
  State->alpha_accept = Proposal.alpha_accept;
}

void lsm_update_Z(LSMModel *Model, LSMState *State)
{
  int R = (State->Z).nrow();
  int C = (State->Z).ncol();
  NumericVector idx = Model->Z_idx - 1;	// R indices start at 1

  LSMState Proposal = copy_lsm_state(State);
  NumericMatrix orig_Z = clone(Proposal.Z);

  for (int i=0; i < idx.size(); ++i) {
    Proposal.Z(idx[i], _) = Proposal.Z(idx[i], _)
      + rnorm(C, 0.0, 1.0/idx.size());
  }

  double new_posterior = (Model->lsm_posterior_fn)(Model, &Proposal);
  double r = R::runif(0, 1);
  double p = exp(new_posterior - Proposal.posterior);

  if (r < p) {
    Proposal.posterior = new_posterior;
    Proposal.Z_accept++;
  }
  else {
    Proposal.Z = orig_Z;
  }

  State->Z = Proposal.Z;
  State->posterior = Proposal.posterior;
  State->Z_accept = Proposal.Z_accept;
}

double log_posterior_logit(LSMModel *Model, LSMState *State)
{
  NumericVector lp = (State->alpha - dist_euclidean(State->Z));
  double post = llik_logit(Model, lp)
    + log_prior_alpha(Model, State)
    + log_prior_Z(Model, State);

  return post;
}

double log_prior_alpha(LSMModel *Model, LSMState *State)
{
  // Matches the prior used by HRH (2002)
  return(R::dgamma(State->alpha, 1.0, 1.0, 1));
}

double log_prior_Z(LSMModel *Model, LSMState *State)
{
  int n = State->Z.ncol();
  NumericVector mu(n);
  NumericMatrix sigma(n);

  mu.fill(0.0);
  sigma.fill_diag(2.0);

  return(sum(dmvnorm(as<arma::mat>(State->Z), as<arma::rowvec>(mu),
		     as<arma::mat>(sigma), true)));
}
