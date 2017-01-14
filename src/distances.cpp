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


#include "dynnet.h"

using namespace Rcpp;

// Euclidean distance

inline double euclidean(NumericVector x, NumericVector y)
{
  double d = 0.0;

  for (int i = 0; i < x.size(); ++i) {
    d += pow(x[i] - y[i], 2.0);
  }

  return(sqrt(d));
}

// [[Rcpp::export(.C_dist_euclidean)]]
NumericVector dist_euclidean(NumericMatrix X)
{
  int nrow = X.nrow();

  NumericVector d(nrow * (nrow-1) / 2);
  int idx = 0;

  for (int j=0; j < (nrow-1); ++j) {
    for (int i=(j+1); i < nrow; ++i) {
      d[idx] = euclidean(X(j,_), X(i,_));
      idx++;
    }
  }

  return(d);
}

// Squared Euclidean distance

inline double euclidean2(NumericVector x, NumericVector y)
{
  double d = 0.0;

  for (int i = 0; i < x.size(); ++i) {
    d += pow(x[i] - y[i], 2.0);
  }

  return(d);
}

// [[Rcpp::export(.C_dist_euclidean2)]]
NumericVector dist_euclidean2(NumericMatrix X)
{
  int nrow = X.nrow();

  NumericVector d(nrow * (nrow-1) / 2);
  int idx = 0;

  for (int j=0; j < (nrow-1); ++j) {
    for (int i=(j+1); i < nrow; ++i) {
      d[idx] = euclidean2(X(j,_), X(i,_));
      idx++;
    }
  }

  return(d);
}

// [[Rcpp::export]]
NumericVector norm_euclidean(NumericMatrix X)
{
  int nrow = X.nrow();
  int ncol = X.ncol();

  NumericVector norm(nrow);

  for (int i = 0; i < nrow; ++i) {
    double sumsqr = 0;
    for (int j = 0; j < ncol; ++j) {
      sumsqr += X(i, j) * X(i, j);
    }
    norm[i] = sqrt(sumsqr);
  }

  return norm;
}

// [[Rcpp::export]]
double distance_penalty(NumericMatrix X)
{
  int nrow = X.nrow();
  int ncol = X.ncol();

  double penalty = 0;

  for (int i = 0; i < nrow; ++i) {
    double sumsqr = 0;
    for (int j = 0; j < ncol; ++j) {
      sumsqr += X(i, j) * X(i, j);
    }
    penalty += log(sqrt(sumsqr));
  }

  return penalty;
}
