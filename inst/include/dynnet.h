#ifndef __DYNNET_H__
#define __DYNNET_H__

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

// structs
struct LSMState
{
    NumericVector beta;
    NumericMatrix Z;
    double posterior;
    int beta_accept;
    int Z_accept;
};

struct LSMModel
{
    NumericVector y;
    NumericMatrix X;		/* model matrix */
    NumericVector Z_idx;	/* indices for Z to be est */
    int k;			/* number of exogenous covariates */
    int d;			/* number of dimensions */
    int burnin;
    int samplesize;
    int interval;
    double (*lsm_posterior_fn)(LSMModel*, LSMState*);
};

// likelihoods
double llik_logit(LSMModel *Model, NumericVector lp);
double llik_poisson(LSMModel *Model, NumericVector lp);

// posteriors
double log_posterior_logit(LSMModel *Model, LSMState *State);

// priors
double log_prior_beta(LSMModel *Model, LSMState *State);
double log_prior_Z(LSMModel *Model, LSMState *State);

// distributions
arma::vec dmvnorm(arma::mat x,
		  arma::rowvec mean,
		  arma::mat sigma,
		  bool logd = false);

arma::mat rmvnorm(int n,
		  arma::vec mu,
		  arma::mat sigma);

// distances
NumericVector dist_euclidean(NumericMatrix X);

// MH
void lsm_update_Z(LSMModel *Model, LSMState *State);
void lsm_update_beta(LSMModel *Model, LSMState *State);
void save_sample(LSMState *State, NumericMatrix *samples, int s);

// message wrappers
void msg_mcmc_iter(int iter, int total, int beta_accept, int Z_accept);

#endif
