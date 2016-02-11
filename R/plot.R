plot.dynnet <- function(network, layout=NULL, ...)
{
    g <- graphs(network)
    t <- length(g)

    if (is.null(layout))
    {
        if (t <= 15) {
            r <- ceiling(t/5)
            c <- ceiling(t/r)
        } else {
            r <- ceiling(15/5)
            c <- ceiling(15/r)
        }
        layout <- c(r, c)
    }

    par(mfrow=layout)
    for (i in seq_along(g))
        plot(g[[i]], ...)
}

plot.dynnetlsm <- function(model, ...)
{
    if (model$method == "MLE") {
        est <- model$estimate$par[-1]            # need a more general way to track the beta coef idx
        est <- matrix(est, ncol=model$k)
    } else if (model$method == "MH") {
        ## do stuff
    }

    G <- model$graph

    all_pos <- insert_ref(est, model$ref, model$k)

    if (model$k == 1)
        all_pos <- cbind(all_pos, 0)

    Y <- as_adj(G)

    plot(G, layout=all_pos, ...)
}
