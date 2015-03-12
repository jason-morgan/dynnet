##' Generate a dynamic network object from a list of adjacency matrices.
##'
##' Generate a dynamic network object from a list of adjacency matrices.
##' @title Dynamic Network from a List of Adjacency Networks
##' @param adj List of adjacency matrixes. Single matrices, representing a
##' static network, are also supported.
##' @param X A \code{data.frame} of nodal covariates. Default is \code{NULL}.
##' @param mode Type of graph. See \code{\link{igraph::graph.adjacency}}.
##' @param weighted Indicates whether the network is weighted
##' (non-binary). Default is \code{NULL}. See
##' \code{\link{igraph::graph.adjacency}} for details.
##' @return A \code{dynnet} object.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
dynnet_adjacency <- function(adj, X=NULL, mode="lower", weighted=NULL)
{
    if (is.matrix(adj) || is.data.frame(adj)) {
        adj <- list(adj)
        periods <- 1
    } else {
        periods <- length(adj)
    }

    adj <- lapply(adj, igraph::graph.adjacency, mode=mode, weighted=weighted)
    nodes <- sapply(adj, igraph::vcount)

    structure(list(nodes=nodes, periods=periods, graphs=adj, X=X,
                   mode=mode, weighted=weighted),
              class="dynnet")
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
