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

plot.lsmfit <- function(model, transform="procrustes", ...)
{
    if (model$method == "MLE") {
        est <- model$estimate$par[-model$beta_idx]
    } else if (model$method == "MH") {
        if (transform == "procrustes" && is.null(model$ref)) {
            est <- colMeans(model$estimate$transformed)
        } else {
            est <- model$estimate$samples[,-model$beta_idx]
            est <- colMeans(est)
        }
    }

    est <- matrix(est, ncol=model$d)

    G <- model$graph

    if (!is.null(model$ref)) {
        all_pos <- insert_ref(est, model$ref, model$d)
    } else {
        all_pos <- est
    }

    if (model$k == 1)
        all_pos <- cbind(all_pos, mean(all_pos))

    plot(G, layout=all_pos, ...)
}

plot_mcmc <- function(model, ...)
{
    if (!is.null(model$estimate$samples))
        plot(model$estimate$samples, ...)
}
