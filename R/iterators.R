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

dmap <- function(n, fn, ...)
{
    output <- vector("list", length=n*(n-1)/2)
    idx <- 1

    for (.J in 1:(n-1)) {
        for (.I in (.J+1):n) {
            environment(fn) <- environment()
            output[[idx]] <- fn(...)
            idx <- idx + 1
        }
    }

    output
}

fordyads <- function(n, fn, ...)
{
    for (.J in 1:(n-1)) {
        for (.I in (.J+1):n) {
            environment(fn) <- environment()
            fn(...)
        }
    }

    invisible()
}
