##' Extract adjacency matrices from dynnet object.
##'
##' Extract adjacency matrices from dynnet object.
##' @title Extract Adjacency Matrices from dynnet Object
##' @param net dynnet network.
##' @param period Numeric vector.
##' @return List of adjacency matrices, unless a single \code{period} is
##' specified, then a matrix is returned.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
get_adjacency <- function(network, period=NULL)
{
    if (is.null(period))
        lapply(network[["graphs"]][], `[`)
    else if (length(period) == 1)
        network[["graphs"]][[period]][]
    else
        network[["graphs"]][period][]
}

get_graph <- function(network, period=NULL)
{
    if (is.null(period))
        network[["graphs"]]
    else if (length(period) == 1)
        network[["graphs"]][[period]]
    else
        network[["graphs"]][period]
}

##' Calculates the number of unique dyads in a network.
##'
##' Calculates the number of unique dyads in a network.
##' @title Calculate the Number of Unique Dyads in a Network
##' @param network A \code{dynnet} network.
##' @return The number of unique dyads in a network.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
ndyads.dynnet <- function(network)
{
    n <- network$nodes * (network$nodes - 1)
    if (!isTRUE(network$directed))
        n <- n / 2

    n
}

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

##' Select a random point from a k-dimensional unit hypersphere.
##'
##' Select a random point from a k-dimensional unit hypersphere.
##' @title Random Point From \code{k}-dimensional Unit Hypersphere
##' @param k Positive integer indicating the dimensions of the hypersphere.
##' @return Numeric vector of coordinates.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
rksphere <- function(k)
{
    x <- rnorm(k)
    (1 / sqrt(x %*% x)) * x
}
