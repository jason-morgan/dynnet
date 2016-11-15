##' @useDynLib dynnet
NULL

##' Generate a dynamic network object from a list of adjacency matrices.
##'
##' Generate a dynamic network object from a list of adjacency matrices.
##' @title Dynamic Network from a List of Adjacency Matrices
##' @param adj List of adjacency matrixes. Single matrices, representing a
##'     static network, are also supported.
##' @param vattr List of data.frames specifying the vertex attributes. Must be
##'     the same length as \code{adj}.
##' @param mode Type of graph. See \code{\link{igraph::graph.adjacency}}.
##' @param weighted Indicates whether the network is weighted
##'     (non-binary). Default is \code{NULL}. See
##'     \code{\link{igraph::graph.adjacency}} for details.
##' @param X A list of \code{data.frame}s recording nodal covariates for each
##'     node and periods. This must either have the same number of elements as
##'     adj or contain a single \code{data.frame}, which will then be used for
##'     each time period. Default is \code{NULL}.
##' @return A \code{dynnet} object.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
dynnet_adjacency <- function(adj, vattr=NULL, mode="lower", weighted=NULL)
{
    if (is.matrix(adj) || is.data.frame(adj)) {
        adj <- list(adj)
        periods <- 1
    } else {
        periods <- length(adj)
    }

    if (!is.null(vattr) && is.data.frame(vattr)) {
        vattr <- list(vattr)
    }

    if (!is.null(vattr) && (length(adj) != length(vattr)))
        stop("Adjacency list and vertex attribute list have different lengths")

    if (!all(sapply(vattr, is.data.frame)))
        stop("Vertex attrbutes must be in the form of a data.frame")

    g     <- lapply(adj, igraph::graph.adjacency, mode=mode, weighted=weighted)
    nodes <- sapply(g, igraph::vcount)

    if (!is.null(vattr))
        g <- mapply(assign_vattr, g=g, vattr=vattr, SIMPLIFY=FALSE)

    structure(list(nodes=nodes, periods=periods, graphs=g,
                   mode=mode, weighted=weighted),
              class="dynnet")
}

assign_vattr <- function(g, vattr)
{
    for (n in colnames(vattr)) {
        vertex_attr(g, name=n) <- vattr[,n]
    }

    g
}

##' Create dynnet Network Object from Another Format
##'
##' Create dynnet Network Object from Another Format
##' @title Create dynnet Network Object from Another Format
##' @param network Network object
##' @param ... Other things
##' @return dynnet network object
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
to_dynnet <- function(network, ...)
{
    UseMethod("to_dynnet")
}

##' @export
##' @rdname to_dynnet
to_dynnet.network <- function(network, ...)
{
    adj <- network::as.matrix.network(network)

    vattr <- network::list.vertex.attributes(network)
    Df <- lapply(vattr, function(a) network::get.vertex.attribute(network, a))
    Df <- as.data.frame(Df)
    colnames(Df) <- vattr

    dynnet_adjacency(adj, vattr=Df)
}

##' The graph object making up the network.
##'
##' The graph object making up the network.
##' @title Extract Graphs from Network
##' @param network dynnet network object.
##' @return List of graphs representing each period in the network.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
graphs <- function(network)
{
    UseMethod("graphs")
}

##' @export
##' @rdname graphs
graphs.dynnet <- function(network)
{
    network$graphs
}

##' The number of periods in the dynamic network.
##'
##' The number of periods in the dynamic network.
##' @title Number of Periods in the Network
##' @param network dynnet network object.
##' @return Number of periods in the (possibly dynamic) network.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
periods <- function(network)
{
    UseMethod("periods")
}

##' @export
##' @rdname periods
periods.dynnet <- function(network)
{
    network$periods
}

##' The number of nodes in the dynamic network.
##'
##' The number of nodes in the dynamic network.
##' @title Number of Nodes in the Network
##' @param network dynnet network object.
##' @return Number of nodes in the (possibly dynamic) network.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
nodes <- function(network)
{
    UseMethod("nodes")
}

##' @export
##' @rdname nodes
nodes.dynnet <- function(network)
{
    network$nodes
}

##' Print dynnet network object
##'
##' Print dynnet network object
##' @title Print dynnet Network Object
##' @param network dynnet network object.
##' @return Prints.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
print <- function(network)
{
    UseMethod("print")
}

##' @export
##' @rdname print
print.dynnet <- function(network)
{
    cat("DYNNET\n")
    cat("  |- Periods:", periods(network), "\n")
    cat("  |- Nodes:  ",   nodes(network), "\n")
}

##' Apply function to each network realization. Currently a wrapper around lapply.
##'
##' Apply function to each network realization.
##' @title Apply Function To Each Network Realization
##' @param network dynnet network object.
##' @param fn Function taking an igraph object.
##' @return List of results of applying fn to each realization of the network.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
gapply <- function(network, fn)
{
    UseMethod("gapply")
}

##' @export
##' @rdname gapply
gapply.dynnet <- function(network, fn)
{
    fn <- match.fun(fn)
    lapply(graphs(network), fn)
}

c.dynnet <- function(..., recursive=FALSE)
{
    nets <- list(...)
    tmp  <- nets[[1]]
    tmp$graphs   <- do.call(c, lapply(nets, graphs))
    tmp$periods  <- length(tmp$graphs)
    tmp$nodes    <- sapply(tmp$graphs, igraph::vcount)
    tmp$mode     <- lapply(nets, `[[`, "mode")
    tmp$weighted <- lapply(nets, `[[`, "weighted")
    tmp
}

##' Get the edge density of each graph in the network
##'
##' Get the edge density of each graph in the network
##' @title Network Edge Density
##' @param network dynnet network object.
##' @return Edge densities
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
edge_density <- function(network)
{
    UseMethod("edge_density")
}

##' @export
##' @rdname edge_density
edge_density.dynnet <- function(network)
{
    gapply(network, igraph::edge_density)
}


procrustes <- function(Z, Zstar)
{
    k <- ncol(Z)
    Z <- scale(Z, center=TRUE, scale=FALSE)
    Ztran <- colMeans(Zstar)
    for (i in 1:k)
        Z[,i] <- Z[,i] + Ztran[i]

    A <- t(Z) %*% (Zstar %*% t(Zstar)) %*% Z
    E <- eigen(A, symmetric=TRUE)
    H <- E$vec %*% diag(sqrt(E$val)) %*% t(E$vec)

    t(t(Zstar) %*% Z %*% solve(H) %*% t(Z))
}
