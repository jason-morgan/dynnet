## ----------------------------------------------------------------------------
## This file is part of dynnet
##
## Copyright (C) 2016 Jason W. Morgan <jason.w.morgan@gmail.com>
##
## dynnet and is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program.  If not, see <http://www.gnu.org/licenses/>.
##
## ----------------------------------------------------------------------------


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
    est <- locations(model, transform=transform)
    G <- model$graph
    plot(G, layout=est, ...)
}

plot_mcmc <- function(model, transformed=FALSE, ...)
{
    if (!is.null(model$estimate$samples)) {
        if (isTRUE(transformed))
            plot(model$estimate$transformed, ...)
        else
            plot(model$estimate$samples, ...)
    }
}

plot_samples <- function(model, nsamp=100, transformed=FALSE, ...)
{
    ## stopifnot(model$method == "MH",
    ##           "this model does not contain any posterior samples")

    ## stopifnot(model$d %in% 1:2,
    ##           "support for models with d>2 coming soon...")

    ## Should this just be another utility function?
    proc_Z <- function(theta, model)
    {
        Z <- decompose_theta(theta, model$beta_idx, model$d)$Z

        if (!is.null(model$ref)) {
            Z <- insert_ref(Z, model$ref, model$d)
        }

        Z
    }

    if (isTRUE(transformed)) {
        idx <- sample(1:nrow(model$estimate$transformed), size=nsamp)
        S <- model$estimate$transformed[idx,]
        Zs <- do.call(rbind, lapply(1:nsamp, function(i) matrix(S[i,], ncol=model$d)))
    } else {
        idx <- sample(1:nrow(model$estimate$samples), size=nsamp)
        S <- model$estimate$samples[idx,]
        Zs <- do.call(rbind, lapply(1:nsamp, function(i) proc_Z(S[i,], model)))
    }

    ## Should add a getter for this.
    n <- vcount(model$graph)

    plot(Zs, col=1:n, ...)
}
