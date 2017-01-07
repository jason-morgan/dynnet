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

// [[Rcpp::export(.C_llik_logit)]]
double C_llik_logit(NumericVector y, NumericVector lp)
{
  double llik = sum(y * lp - log1p(exp(lp)));
  return(llik);
}

double llik_logit(LSMModel *Model, NumericVector lp)
{
  double llik = sum(Model->y * lp - log1p(exp(lp)));

  // Left for reference. The above is faster, but this makes it clear what is
  // being computed.
  // for (int i = 0; i < n; ++i) {
  //   double p = R::plogis(lp[i], 0.0, 1.0, 1.0, false);
  //   llik += R::dbinom(Model->y[i], 1.0, p, true);
  // }

  return(llik);
}

// [[Rcpp::export(.C_llik_poisson)]]
double C_llik_poisson(NumericVector y, NumericVector lp)
{
  int n = y.size();
  double llik = 0;

  for (int i = 0; i < n; ++i) {
    double lambda = exp(lp[i]);
    llik += R::dpois(y[i], lambda, true);
  }

  return(llik);
}

double llik_poisson(LSMModel *Model, NumericVector lp)
{
  int n = Model->y.size();
  double llik = 0;

  for (int i = 0; i < n; ++i) {
    double lambda = exp(lp[i]);
    llik += R::dpois(Model->y[i], lambda, true);
  }

  return(llik);
}
