plot.dynnet <- function(network, ...)
{
    g <- graphs(network)
    t <- length(g)

    if (t <= 15) {
        r <- ceiling(t/5)
        c <- ceiling(t/r)
    } else {
        r <- ceiling(15/5)
        c <- ceiling(15/r)
    }

    par(mfrow=c(r, c))
    for (i in seq_along(g))
        plot(g[[i]], ...)
}
