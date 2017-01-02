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

// simple adaptive proposals
void adapt_Z_proposal_sd(LSMState *State, int niter, int Z_accept_last)
{
  double rate = (State->Z_accept - Z_accept_last) / 2000.0;

  if ((rate > 0.30) && (State->Z_proposal_sd < 10.0)) {
    State->Z_proposal_sd = State->Z_proposal_sd * 1.20;
    Rcpp::Rcout << "     => Increasing Z proposal SD to "
		<< State->Z_proposal_sd
		<< " ( last Z acceptance rate: "
		<< rate
		<< " )"
		<< std::endl;
  }

  if ((rate < 0.20) && (State->Z_proposal_sd > 0.001)) {
    State->Z_proposal_sd = State->Z_proposal_sd * 0.80;
    Rcpp::Rcout << "     => Decreasing Z proposal SD to "
		<< State->Z_proposal_sd
		<< " ( last Z acceptance rate: "
		<< rate
		<< " )"
		<< std::endl;

  }
}
