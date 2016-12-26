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


##' Simple linear trajectory.
##'
##' Simple linear trajectory.
##' @title Linear Trajectory
##' @param period Integer specifying the time period for which to calculate
##' the latent positions.
##' @param pos Matrix of latent locations.
##' @param traj Matrix of trajectories.
##' @param ... Additional parameters to pass. Not used in this trajectory
##' function.
##' @return Matrix of latent positions.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
trajectory_linear <- function(period, pos, traj, ...)
{
    pos + (period - 1) * traj
}
