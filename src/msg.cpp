#include "dynnet.h"

using namespace Rcpp;

// Wrapper around a message printer
void msg_mcmc_iter(int iter, int total,
		   int beta_accept, int Z_accept)
{
  Rcpp::Rcout << "iter => "
	      << iter
	      << " (" << (iter * 100) / total << "%)  "
	      << "( beta accept: "
	      << (double(beta_accept) / double(iter)) * 100
	      << "%  Z accept: "
	      << (double(Z_accept) / double(iter)) * 100
	      << "% )"
	      << std::endl;
}
