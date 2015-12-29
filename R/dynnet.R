##' Generate a dynamic network object from a list of adjacency matrices.
##'
##' Generate a dynamic network object from a list of adjacency matrices.
##' @title Dynamic Network from a List of Adjacency Matrices
##' @param adj List of adjacency matrixes. Single matrices, representing a
##' static network, are also supported.
##' @param X A list of \code{data.frame}s recording nodal covariates for each
##' node and periods. This must either have the same number of elements as adj
##' or contain a single \code{data.frame}, which will then be used for each time
##' period. Default is \code{NULL}.
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

    adj   <- lapply(adj, igraph::graph.adjacency, mode=mode, weighted=weighted)
    nodes <- sapply(adj, igraph::vcount)

    structure(list(nodes=nodes, periods=periods, graphs=adj, X=X,
                   mode=mode, weighted=weighted),
              class="dynnet")
}

##' The number of periods in the dynamic network.
##'
##' The number of periods in the dynamic network.
##' @title Number of Periods in the Network
##' @param net dynnet network object.
##' @return Number of periods in the (possibly dynamic) network.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
periods.dynnet <- function(net)
{
    net$periods
}
