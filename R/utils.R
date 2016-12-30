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
##' @keywords internal
list_all_equal <- function(lst)
{
    length(unique(lst)) == 1
}

##' Decompose estimated theta vector into its beta and Z components.
##'
##' This function is meant as the go-to for decomposing theta (such as a single
##' sample from the posterior or the vector of estimates from optim). The list
##' returned includes a vector containing the coefficients on the exogenous
##' covariates and a positions matrix, Z.
##' @title Decompose \code{theta} into the
##' @param theta Vector to decompose.
##' @param beta_idx Integer vector indicating the
##' @param d Dimensions of the latent space.
##' @return List of beta, Z
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @keywords internal
decompose_theta <- function(theta, beta_idx, d)
{
    list(beta=theta[beta_idx],
         Z=matrix(theta[-beta_idx], ncol=d))
}

##' Insert the reference vertex locations into the latent positions matrix.
##'
##' This calls an underlying \code{C} function, \code{.C_insert_ref}, since it's
##' critical to the speed of the MLE implementation.
##' @title Insert Reference Locations into Latent Positions Matrix
##' @param pos \code{n-r} by \code{d} positions matrix, where r is the number of
##' @param ref Reference object as specified in an latent space model
##' @param d Dimensions of the latent space
##' @return Latent positions matrix with reference locations inserted.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @keywords internal
insert_ref <- function(pos, ref, d)
{
    n <- nrow(pos) + length(ref$idx)
    .C_insert_ref(ref$idx, ref$pos, (1:n)[-ref$idx], pos)
}
