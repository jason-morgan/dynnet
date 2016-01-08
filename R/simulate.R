##' Simulate dynammic and static networks.
##'
##' Simulate dynammic and static networks.
##' @title Simulate Dynamic and Static Networks
##' @param n Positive integer. Number of nodes in the network.
##' @param k Positive integer. Number of dimensions in the latent space.
##' @param periods Positive integer. Number of periods in the dynnamic
##'     network. If \code{k=1}, a static network is generated.
##' @param vattr List of vertex attributes. One data.frame per period.
##' @param ref_fn Function to generate the reference positions. Function must
##'     return a list reference positions and indices.
##' @param seed Set the seed before generating the networks to assure
##'     replicability. Default is \code{NULL}.
##' @param family Character string in \code{c("bernoulli", "poisson")}.
##' @param ... Further parameters to be passed to subsequent functions.
##' @return A \code{dynsim} object.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
simulate_dynnet <- function(n, k=1, periods=1,
                             vattr   = NULL,
                             ref_fn  = default_ref_pos(k),
                             seed=NULL, family="bernoulli", ...)
{
    if (!is.null(seed))
        set.seed(seed)

    if (n < k+1)
        stop("Number of nodes cannot be less than the number of reference",
             " nodes (k+1)")

    if (!is.null(vattr) & length(vattr) != periods)
        stop("vattr must has the same length as the number of periods to be",
             " simulated.")

    ## Generate starting position matrix
    N   <- split_n(n-nrow(ref$pos), g)
    pos <- generate_latent_values(N, ref$pos, groups$mean, groups$sigma)

    ## Generate trajectory matrix
    if (periods > 1) {
        ## Generate trajectory matrix
        N    <- split_n(n-nrow(ref$traj), g)
        traj <- generate_latent_values(N, ref$traj, groups$traj,
                                       groups$sigma_traj)
        mv   <- calc_positions(1:periods, pos, traj, fn=traj_fn, ...)
    } else {
        traj <- NULL
        mv   <- list(pos)
    }

    ## Generate adjacency matrices.
    Xb  <- mapply(`%*%`, X, list(beta), SIMPLIFY=FALSE)
    adj <- mapply(generate_adjacency, mv, Xb, family=family, SIMPLIFY=FALSE)

    ## Create dynnet object.
    wght <- if (family == "bernoulli") NULL else TRUE
    SIM <- dynnet_adjacency(adj, weighted=wght)

    structure(list(DGP=list(positions=pos, trajectories=traj,
                            ref=list(pos=ref$pos, idx=1:(k+1)),
                            beta=beta, family=family, seed=seed),
                   SIM=SIM),
              class="dynnetsim")
}

##' Generate latent positions or trajectories for nodes in the network.
##'
##' Generate latent positions or trajectories for nodes in the network.
##' @title Generate Latent Positions or Trajectories
##' @param n List of positive integers indicating the number of nodes in each
##' cluster, including the reference unit.
##' @param ref List of vectors indicating the location of the reference units.
##' @param mean List of mean vectors indicating the mean for each cluster of
##' nodes.
##' @param sigma List of covariance matrices, one for each cluster of nodes.
##' @return A matrix of positions or trajectories in the latent space.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
generate_latent_values <- function(n, ref, mean, sigma)
{
    mean <- split(mean, 1:nrow(mean))

    x <- mapply(generate_group, n=n, mean=mean, sigma=sigma, SIMPLIFY=FALSE)
    x <- rbind(ref, do.call(rbind, x))
    colnames(x) <- paste0("d", 1:ncol(x))
    rownames(x) <- c(paste0("r", 1:nrow(ref)),
                     paste0("n", 1:(nrow(x)-nrow(ref))))
    x
}

##' Produces a function defining the data generating process for ties in the
##' network.
##'
##' Produces a function defining the data generating process for ties in the
##' network.
##' @title Network Data Generating Process
##' @param family String in \code{c("logit", "poisson")} specifying the data
##' generating process of ties in the network.
##' @return A function taking single argument, a linear predictor.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
dgp_fn <- function(family)
{
    switch(family,
           "bernoulli" = function(lp) { rbinom(length(lp), 1, plogis(lp)) },
           "poisson"   = function(lp) { rpois(length(lp), lambda=exp(lp)) },
           stop("unknown family"))
}

##' Generate an adjacency matrix from a matrix of latent positions, nodal
##' covariates, and a vector of coefficients on those covariates.
##'
##' Generate an adjacency matrix from a matrix of latent positions, nodal
##' covariates, and a vector of coefficients on those covariates.
##' @title Generate Adjacency Matrix from a Set of Latent Positions, Covariate
##' Profile, and Coefficient Vector
##' @param pos Matrix of latent positions.
##' @param Xb Vector of linear predictors for each dyad.
##' @param family String indicating the data generating process to follow.
##' @param ... Additional parameters to be passed to the distance function.
##' @return An adjacency matrix. These matrices are of a type produced by the
##' \code{\link{Matrix}} package.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
generate_adjacency <- function(pos, Xb, family=NULL, ...)
{
    A <- matrix(0, ncol=nrow(pos), nrow=nrow(pos))
    d <- as_distance_vector(pos, ...)
    A[lower.tri(A)] <- dgp_fn(family)(Xb - d)
    Matrix::Matrix(A)
}

##' Divides an integer into \code{g} approximately equal integers and returns as
##' a list.
##'
##' Divides an integer into \code{g} approximately equal integers and returns as
##' a list.
##' @title Splits \code{n} into \code{g} Groups
##' @param n Integer indicating the total number of nodes in the network.
##' @param g Integer indicating the number of groups.
##' @return List of integers of approximately the same size.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
split_n <- function(n, g)
{
    N <- rep(floor(n / g), g)
    r <- n %% g
    if (r > 0) N[1:r] <- N[1:r] + 1
    split(N, 1:g)
}

##' Create random default reference node locations on a \code{k}-dimensional
##' unit hypersphere.
##'
##' Create random default reference node locations on a \code{k}-dimensional
##' unit hypersphere. The minimum distance between any two nodes is guaranteed
##' to be greater than \code{sqrt(2) / scale} where \code{scale = k-1} when
##' \code{k > 2}. In 2 dimensions, nodes will never be located within the same
##' quadrant. \code{scale} is used for larger dimensions to greatly reduce
##' reference position generation time (using \code{sqrt(2)} meant that many
##' sets of reference positions were being rejected).
##' @title Generate Default Reference Locations
##' @param k Positive integer indicating the dimensions of the latent social
##' space.
##' @return An \code{(k+1)}-by-\code{k} numeric matrix of latent locations on
##' the \code{k}-dimensional unit hypersphere.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
default_ref_pos <- function(k)
{
    scale <- if (k < 3) 1 else k-1

    if (k == 1) {
        x <- matrix(c(-1, 1), ncol=1)
    } else {
        x <- c(0, 0)
        while (min(dist(x)) < sqrt(2)/scale) {
            x <- t(replicate(k+1, rksphere(k)))
        }
    }
    x
}

##' Extract latent positions from simulated dynamic network.
##'
##' Extract latent positions from simulated dynamic network.
##' @title Extract Latent Positions from Simulated Network
##' @param sim Simulated network of class \code{\link{dynnetsim}}.
##' @return Matrix of latent positions.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
positions <- function(sim) { sim$DGP$positions }

##' Extract latent trajectories from simulated dynamic network.
##'
##' Extract latent trajectories from simulated dynamic network.
##' @title Extract Latent Trajectories from Simulated Network
##' @param sim Simulated network of class \code{\link{dynnetsim}}.
##' @return Matrix of latent trajectories.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
trajectories <- function(sim) { sim$DGP$trajectories }
