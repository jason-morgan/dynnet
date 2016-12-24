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

plot.lsmfit <- function(model, ...)
{
    if (model$method == "MLE") {
        est <- model$estimate$par[-model$beta_idx]
    } else if (model$method == "MH") {
        est <- model$estimate$samples[,-model$beta_idx]
        est <- colMeans(est)
    }

    est <- matrix(est, ncol=model$d)

    G <- model$graph

    if (!is.null(model$ref))
        all_pos <- insert_ref(est, model$ref, model$d)
    else
        all_pos <- est

    if (model$k == 1)
        all_pos <- cbind(all_pos, 0)

    plot(G, layout=all_pos, ...)
}
