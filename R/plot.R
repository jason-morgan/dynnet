## plot <- function(network, layout=NULL, ...)
## {
##     UseMethod("plot")
## }

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

plot.dynnetlsm <- function(model)
{
    if (model$method == "MLE") {
        est <- model$estimate$par[-1]            # need a more general way to track the beta coef idx
        est <- matrix(est, ncol=model$k)
    } else if (model$method == "MH") {
        ## do stuff
    }

    all_pos <- insert_ref(est, model$ref, model$k)
    Y <- as_adj(model$graph)

    plot(all_pos, pch=19, col="white")

    plot_line <- function(Y, pos)
    {
        if (Y[.I,.J] != 0)
            lines(pos[c(.I,.J),1], pos[c(.I,.J),2])
    }

    fordyads(nrow(all_pos), plot_line, Y=Y, pos=all_pos)

    points(est, pch=19)
    points(model$ref$pos, col="red", cex=1.25)
    points(model$ref$pos, col="red", pch=19, cex=0.5)
}
