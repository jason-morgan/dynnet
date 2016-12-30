## ----------------------------------------------------------------------------
## This file is part of dynnet
##
## Copyright (C) 2016 Jason W. Morgan <jason.w.morgan@gmail.com>
##
## dynnet and is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program.  If not, see <http://www.gnu.org/licenses/>.
##
## ----------------------------------------------------------------------------


##' Converts a matrix of latent locations to a vector of pairwise distances.
##'
##' Converts a matrix of latent locations for each of the \code{n} nodes to a
##' column vector of pairwise distances.
##' @title Convert Latent Locations to a Vector of Distances
##' @param Z Matrix representing latent locations.
##' @param ... Additional parameters to be passed to \code{\link{dist}}.
##' @return Column vector of distances as contained in the \code{lower.tri(d)},
##' where \code{d} is the distance matrix.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
as_distance_vector <- function(Z, ...)
{
    matrix(dist(Z, ...))
}

##' Calculate node positions for a specific time period.
##'
##' Calculate node positions for a specific time period.
##' @title Calculate Node Postitions
##' @param periods Integer vector of time periods for which to calculate the
##' position.
##' @param pos Matrix of latent locations.
##' @param traj Matrix of trajectories.
##' @param fn Function taking a time period, latent starting positions, and
##' trajectories, which returns a matrix of latent position. Default is the
##' simple linear trajectory model.
##' @param ... List of parameters to pass to the trajectory function.
##' @return List of latent locations, one for each time period.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
calc_positions <- function(periods, pos, traj, fn=trajectory_linear, ...)
{
    lapply(periods, fn, pos, traj, ...)
}

##' Checks whether all elements of the supplied list are equal.
##'
##' Checks whether all elements of the supplied list are equal.
##' @title Check if All List Elements are Equal
##' @param lst List.
##' @return Logical indicating whether or not all values of the list are equal.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
list_all_equal <- function(lst)
{
    length(unique(lst)) == 1
}

insert_ref <- function(pos, ref, d)
{
    n <- nrow(pos) + length(ref$idx)
    .C_insert_ref(ref$idx, ref$pos, (1:n)[-ref$idx], pos)
}
