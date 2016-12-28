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
