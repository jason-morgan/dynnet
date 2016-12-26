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


vcov.lsmfit <- function(object, ...)
{
    solve(-object$estimate$hessian)
}

coef.lsmfit <- function(object, ...)
{
    idx <- object$beta_idx
    if (object$method == "MLE") {
        est <- object$estimate$par[idx]
    } else if (object$method == "MH") {
        if (length(idx) == 1)
            est <- mean(object$estimate$samples[, idx])
        else
            est <- colMeans(object$estimate$samples[, idx])
    }

    names(est) <- object$beta_names
    est
}

summary.lsmfit <- function(object, ...)
{
    idx <- object$beta_idx
    est <- coef(object)

    if (object$method == "MLE") {
        se  <- sqrt(diag(vcov(object))[idx])
        tbl <- data.frame("Estimate"=est, "SE"=se)
    } else if (object$method == "MH") {
        if (length(idx) == 1)
            se <- sd(object$estimate$samples[, idx])
        else
            se <- apply(object$estimate$samples[, idx], 2, sd)

        tbl <- data.frame("Posterior Mean"=est, "SD"=se)
    }

    rownames(tbl) <- object$beta_names
    tbl
}

logLik.lsmfit <- function(object, ...)
{
    if (object$method == "MLE") {
        llik <- object$estimate$value
    } else if (object$method == "MH") {
        llik <- NULL
    }

    llik
}
