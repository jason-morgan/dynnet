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


absdiff <- function(graph, vattr)
{
    x <- vertex_attr(graph, vattr)
    if (is.factor(x))
        stop("variable ", vattr,
             " is a factor. absdiff supports real values only")

    a <- outer(x, x, function(x1, x2) abs(x1 - x2))
    a[lower.tri(a)]
}

covmatch <- function(graph, vattr)
{
    x <- vertex_attr(graph, vattr)
    if (!is.factor(x) || is.character(x)) {
        warning("variable ", vattr,
                " is neither a factor or character vector. ",
                " Coercing to factor (you may not want this).",
                immediate.=TRUE)
        x <- factor(x)
    }

    a <- outer(x, x, function(x1, x2) x1 == x2)
    a[lower.tri(a)]
}
