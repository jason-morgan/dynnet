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


## This file contains utility functions that mirror some of those found in the
## igraph package.

vcount <- function(object, ...)
{
    UseMethod("vcount")
}

##' @export
##' @rdname vcount
vcount.igraph <- function(object, ...)
{
    igraph::vcount(object)
}

##' @export
##' @rdname vcount
vcount.dynnet <- function(object, ...)
{
    gapply(object, igraph::vcount)
}
