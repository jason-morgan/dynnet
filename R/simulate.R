##' Simulate dynammic and static networks.
##'
##' Simulate dynammic and static networks.
##' @title Simulate Dynamic and Static Networks
##' @param n Positive integer. Number of nodes in the network.
##' @param k Positive integer. Number of dimensions in the latent space.
##' @param periods Positive integer. Number of periods in the dynnamic
##' network. If \code{k=1}, a static network is generated.
##' @param groups Positive integer. Number of groups, or clusters, of nodes in
##' the network.
##' @param X List of model matrices of covariates of size \code{N}-by-\code{m},
##' where \code{N} is the \emph{number of dyads} in the network and \code{m} is
##' the number of covariates. If \code{NULL} an intercept is used. If a single
##' matrix is provided, this model matrix will be used for all periods.
##' @param beta Vector of coefficients of length \code{m}. Default is
##' \code{c(0.5)} for the coefficient on the intercept.
##' @param ref_pos Matrix representing the positions of the reference units in
##' the latent space. Should be of dimensions \code{k+1}-by-\code{k}. Default is
##' to randomly choose \code{k+1} points on a \code{k}-dimensional hypersphere.
##' @param mean_pos Matrix representing the mean position of the groups. Should
##' be of dimensions \code{g}-by-\code{k}. Default is to use the origin for all
##' groups.
##' @param sigma_pos List of \code{g} covariance matrices used to generate the
##' locations of nodes in each group. Default is a list of \code{g}
##' \code{k}-by-\code{k} identity matrices.
##' @param ref_traj Matrix representing the trajectories of the reference units
##' through the latent space. Should be of dimensions
##' \code{k+1}-by-\code{k}. For identification, the locations of the reference
##' units should ideally be fixed.
##' @param mean_traj Matrix representing the mean of the node
##' trajectories. Default is a \code{k}-by-\code{g} matrix of zeros.
##' @param sigma_traj List of \code{g} covariance matrices used to generate the
##' trajectories. Default is a list of \code{g} \code{k}-by-\code{k} identity
##' matrices. Note: this may imply too much movement in the nodes.
##' @param traj_fn Function taking a positive integer representing the period
##' for which to generate the network as well as a position and trajectory
##' matrix. Default is \code{\link{linear_trajectory}}.
##' @param seed Set the seed before generating the networks to assure
##' replicability. Default is \code{NULL}.
##' @param family Character string in \code{c("logit", "poisson")}.
##' @param ... Further parameters to be passed to subsequent functions.
##' @return A \code{dynsim} object.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
simulate_network <- function(n, k=1, g=1, periods=1,
                             X       = NULL,
                             beta    = 0.5,
                             ref     = list(pos=default_ref_pos(k),
                                            idx=1:(k+1),
                                            traj=matrix(0, ncol=k, nrow=k+1),
                                            sigma_traj=NULL),
                             groups  = list(mean=matrix(0, ncol=k, nrow=g),
                                            sigma=rep(list(diag(k)), g),
                                            traj=matrix(0, ncol=k, nrow=g),
                                            sigma_traj = rep(list(diag(k)),
                                                             g)),
                             traj_fn = trajectory_linear,
                             seed=NULL, family="logit", ...)
{
    if (is.null(X)) {
        X <- rep(list(matrix(1, n*(n-1)/2)), periods)
    } else if (is.matrix(X)) {
        X <- rep(list(X), periods)
    }

    if (!is.null(seed))
        set.seed(seed)

    if (n < k+1)
        stop("Number of nodes cannot be less than the number of reference",
             " nodes (k+1)")

    if (!all(dim(ref$pos) == c(k+1, k)))
        stop("Reference position matrix (ref$pos) has the wrong dimensions.",
             " Should be of\n  dimensions (k+1, k).")

    if (!all(dim(groups$mean) == c(g, k)))
        stop("Group position matrix (groups$mean) has the wrong dimensions.",
             " should be of\n  dimensions (g, k).")

    if (g != length(groups$sigma))
        stop("Wrong number of group covariances supplied. Length of",
             " groups$sigma should be equal\n  to the number of groups.")

    if (list_all_equal(split(ref$pos, 1:nrow(ref$pos))))
        warning("Reference positions all equal. You probably don't want this.")

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
    wght <- if (family == "logit") NULL else TRUE
    SIM <- dynnet_adjacency(adj, X=X, weighted=wght)

    structure(list(DGP=list(positions=pos, trajectories=traj,
                            ref=list(pos=ref$pos, idx=1:(k+1)),
                            beta=beta, family=family, seed=seed),
                   SIM=SIM),
              class="dynnetsim")
}

##' Generate a cluster of nodes in \code{k}-dimensional latent space.
##'
##' Generate a cluster of nodes in \code{k}-dimensional latent space.
##' @title Generate a Cluster of Nodes
##' @param n Positive integer. The number of nodes in the cluster.
##' @param mean Mean vector. Default is the a zero vector of length \code{k}.
##' @param sigma Covariance matrix. Default is a \code{k}-by-\code{k} identity
##' matrix.
##' @return Numeric matrix with dimensions of \code{(n-1)}-by-\code{k}.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
generate_group <- function(n, mean, sigma)
{
    mvtnorm::rmvnorm(n, mean=mean, sigma=sigma)
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
           "logit"   = function(lp) { rbinom(length(lp), 1, plogis(lp)) },
           "poisson" = function(lp) { rpois(length(lp), lambda=exp(lp)) },
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
