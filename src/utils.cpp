#include "dynnet.h"

using namespace Rcpp;

// [[Rcpp::export(.C_insert_ref)]]
NumericMatrix C_insert_ref(IntegerVector ref_idx,
			   NumericMatrix ref_pos,
			   IntegerVector est_idx,
			   NumericMatrix est_pos)
{
  NumericMatrix result(ref_idx.size() + est_idx.size(),
		       ref_pos.ncol());

  // R indices start at 1
  ref_idx = ref_idx - 1;
  est_idx = est_idx - 1;

  for (int i = 0; i < ref_idx.size(); ++i) {
    result(ref_idx[i],_) = ref_pos(i,_);
  }

  for (int i = 0; i < est_idx.size(); ++i) {
    result(est_idx[i],_) = est_pos(i,_);
  }

  return(result);
}
