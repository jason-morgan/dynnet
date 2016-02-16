#ifndef __DYNNET_H__
#define __DYNNET_H__

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

// structs
struct LSMModel
{
    NumericVector y;
    NumericMatrix X;		/* model matrix */
    NumericVector Z_idx;	/* indices for Z to be est */
    int k;			/* number of dimensions */
    int burnin;
    int samplesize;
    int interval;
};

struct LSMState
{
    double alpha;
    NumericVector beta;
    NumericMatrix Z;
    double posterior;
    int alpha_accept;
    int beta_accept;
    int Z_accept;
};

// likelihoods
double C_llik_logit(NumericVector y, NumericVector lp);
double C_llik_poisson(NumericVector y, NumericVector lp);

// posteriors
double C_lsm_posterior(LSMModel *Model, LSMState *State);
double C_log_posterior_logit(NumericVector y, NumericVector lp,
			     double alpha, NumericMatrix Z);

// distributions
arma::vec C_dmvnorm(arma::mat x,
		    arma::rowvec mean,
		    arma::mat sigma,
		    bool logd = false);

// distances
NumericVector C_dist_euclidean(NumericMatrix X);

// MH
void C_lsm_update_Z(LSMModel *Model, LSMState *State);
void C_lsm_update_alpha(LSMModel *Model, LSMState *State);
void save_sample(LSMState *State, NumericMatrix *samples, int s);

#endif
