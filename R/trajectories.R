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
