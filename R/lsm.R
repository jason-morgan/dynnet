lsm <- function(network, k=1, period=1, ref=NULL, family="logit", method="MLE",
                fit=TRUE, verbose=TRUE)
{
    ## Construct ref
    ## if (is.null(ref)) {
    ##     ref <- list(idx=1:(k+1), pos=default_ref_pos(k))
    ## }
    if (is.null(ref))
        stop("lsm: reference node information not provided", call.=FALSE)

    Y   <- get_graph(network, period=period)
    deg <- igraph::degree(Y)
    Y   <- igraph::delete.vertices(Y, deg == 0)

    if (sum(deg == 0) > 0)
        warning("lsm: isolate removed from network", call.=FALSE)

    model <- list(network=network, period=1, k=k, ref=ref, family=family,
                  dropped=which(deg == 0))

    ## Set covariate matrix
    ndyads <- nrow(Y[]) * (nrow(Y[]) - 1) / 2
    model$X <- matrix(1, nrow=ndyads, ncol=1)

    ## Starting values
    model$start <- init_coef(Y, model$X, ref$pos, ref$idx)

    if (isTRUE(verbose))
        cat("Starting values:", model$start, "\n")

    if (isTRUE(fit))
        lsm_MLE(model, verbose=verbose)
    else
        model
}
