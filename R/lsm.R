lsm <- function(network, k=1, period=1, ref=NULL, family="logit", method="MLE",
                fit=TRUE, verbose=TRUE)
{
    if (is.null(ref))
        stop("lsm: reference node information not provided", call.=FALSE)

    Y   <- get_graph(network, period=period)
    deg <- igraph::degree(Y)
    Y   <- igraph::delete.vertices(Y, deg == 0)

    if (sum(deg == 0) > 0)
        warning("lsm: isolate removed from network", call.=FALSE)

    model <- list(network=network, period=1, k=k, ref=ref, family=family,
                  dropped=which(deg == 0))

    model
}
