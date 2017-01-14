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
    NumericVector old_beta;
    NumericMatrix Z;
    NumericMatrix old_Z;
    NumericVector Xb;
    NumericVector old_Xb;
    NumericVector dist;
    NumericVector old_dist;
    double lpr_graph;
    double posterior;
    double beta_log_prior;
    int beta_accept;
    double Z_log_prior;
    int Z_accept;
    double Z_proposal_sd;
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
    double (*lsm_lpr_fn)(LSMModel*, NumericVector lp);
    NumericVector (*lsm_dist_fn)(NumericMatrix X);
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
NumericVector dist_euclidean2(NumericMatrix X);

// MH
void lsm_update_Z(LSMModel *Model, LSMState *State);
void lsm_update_beta(LSMModel *Model, LSMState *State);
void save_sample(LSMState *State, NumericMatrix *samples, int s);

// Adaptive proposals
void adapt_Z_proposal_sd(LSMState *State, int niter, int Z_accept_last);

// message wrappers
void msg_mcmc_iter(int iter, int total, int beta_accept, int Z_accept);

#endif
