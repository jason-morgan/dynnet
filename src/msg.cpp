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
