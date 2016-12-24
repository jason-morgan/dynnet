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
