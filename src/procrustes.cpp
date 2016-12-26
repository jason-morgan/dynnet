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

// #include <RcppArmadillo.h>

// // [[Rcpp::depends(RcppArmadillo)]]

// using namespace Rcpp;
// using namespace arma;

// // [[Rcpp::export(.C_procrustes)]]
// arma::mat procrustes(arma::mat Z, arma::mat Zstar)
// {
//   int k = Z.n_cols;
//   arma::mat Z_mu = mean(Z,0);
//   arma::mat Zstar_mu = mean(Zstar,0);

//   for (int j=0; j < k; ++j) {
//     Z.col(j) = Z.col(j) - Z_mu(0,j) + Zstar_mu(0,j);
//   }

//   arma::mat A = Z.t() * (Zstar * Zstar.t()) * Z;

//   arma::vec eigval;
//   arma::mat eigvec;
//   arma::eig_sym(eigval, eigvec, A);

//   arma::mat H = eigvec %*% diag(sqrt(eigval)) %*% eigvec.t();

//   t(t(Zstar) %*% Z %*% solve(H) %*% t(Z))

//   return(Z);
// }
