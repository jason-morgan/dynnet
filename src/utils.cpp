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

// [[Rcpp::export(.C_insert_ref)]]
NumericMatrix insert_ref(IntegerVector ref_idx,
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
