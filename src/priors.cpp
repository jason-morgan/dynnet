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
#include "dynnet.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


double log_prior_beta(LSMModel *Model, LSMState *State)
{
  double p = 0.0;

  if (Model->k == 1) {
    // Matches the prior used by HRH (2002)
    p = R::dgamma((State->beta)[0], 1.0, 1.0, 1);
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
  sigma.fill_diag(5.0);

  return(sum(dmvnorm(as<arma::mat>(State->Z),
		     as<arma::rowvec>(mu),
		     as<arma::mat>(sigma), true)));
}
