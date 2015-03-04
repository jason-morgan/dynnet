##' Generate a dynamic network object from a list of adjacency matrices.
##'
##' Generate a dynamic network object from a list of adjacency matrices.
##' @title Dynamic Network from a List of Adjacency Networks
##' @param adj List of adjacency matrixes. Single matrices, representing a
##' static network, are also supported.
##' @param X A \code{data.frame} of nodal covariates. Default is \code{NULL}.
##' @param weighted Indicates whether the network is weighted (non-binary). If a
##' weighted network is provided, but \code{weighted} is not set to \code{TRUE},
##' all elements of the matrices in \code{adj} will be converted to binary
##' values, with any value greater than \code{0} equal to \code{1} and all other
##' set to \code{0}. Default is \code{FALSE}.
##' @param directed Indicates whether the network is directed. Default is
##' \code{FALSE}.
##' @param bipartite Indicates whether the network is bipartite. Default is
##' \code{FALSE}.
##' @return A \code{dynnet} object.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
dynnet_adjacency <- function(adj, X=NULL, weighted=FALSE, directed=FALSE,
                             bipartite=FALSE)
{
    if (is.matrix(adj)) {
        adj <- list(adj)
        periods <- 1
    } else {
        periods <- length(adj)
    }

    if (isTRUE(bipartite))
        nodes <- c(nrow(adj[[1]]), ncol(adj[[1]]))
    else
        nodes <- nrow(adj[[1]])

    ## Force all entries to be binary if weighted is FALSE
    if (!isTRUE(weighted))
        adj <- lapply(adj,
                      function(m)
                      {
                          m[m > 0]  <- 1
                          m[m <= 0] <- 0
                          m
                      })

    structure(list(nodes=nodes, periods=periods, adj=adj, X=X,
                   weighted=weighted, directed=directed, bipartite=bipartite),
              class="dynnet")
}

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
get_adjacency <- function(net, period=NULL)
{
    if (is.null(period))
        net[["adj"]]
    else if (length(period) == 1)
        net[["adj"]][[period]]
    else
        net[["adj"]][period]
}


## create_stan_data <- function()
## {
##     ## int<lower=1> n;            // number of nodes
##     ## int<lower=1> N;            // number of dyads
##     ## int<lower=1> K;            // number of dimensions in latent space
##     ## int<lower=0> y[N];         // outcome
##     ## real x[N];                 // one covariate

##     dta <- data.frame(t_idx=numeric(), node1_idx=numeric(), node2_idx=numeric(),
##                       D=numeric())

##     for (t in 1:T) {
##         d <- as.matrix(dist(pos + (t-1) * traj))
##         node1_idx <- col(d)[lower.tri(d)]
##         node2_idx <- row(d)[lower.tri(d)]
##         dta <- rbind(dta, data.frame(t_idx=t, node1_idx, node2_idx,
##                                      D=d[lower.tri(d)]))
##     }

## }
