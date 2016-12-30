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

##' Extract the estimated latent locations for a latent space model.
##'
##' This function returns the estimated latent locations for a latent space
##' model.
##' @title Estimated Locations from Latent Space Model
##' @param object A \code{\link{lsmfit}} object.
##' @param ... Addition parameters to pass.
##' @return A matrix of estimated locations. For a model fit via MLE, this is
##'     the MLE estimates. For models with with the Metropolis-Hastings
##'     algorithm, the posterior means for the locations (possibly transformed
##'     via a Procrustes step, if that was the method of identification
##'     specified.)
##' @author Jason Morgan
locations <- function(object, ...)
{
    UseMethod("locations")
}

##' @export
##' @rdname locations
locations.lsmfit <- function(object, transform="procrustes", ...)
{
    if (object$method == "MLE") {
        est <- object$estimate$par[-object$beta_idx]
    } else if (object$method == "MH") {
        if (transform == "procrustes" && is.null(object$ref)) {
            est <- colMeans(object$estimate$transformed)
        } else {
            est <- object$estimate$samples[,-object$beta_idx]
            est <- colMeans(est)
        }
    }

    est <- matrix(est, ncol=object$d)

    if (!is.null(object$ref))
        est <- insert_ref(est, object$ref, object$d)

    if (object$d == 1)
        est <- cbind(est, mean(all_pos))

    rownames(est) <- vertex_attr(object$graph, "name")
    est
}

## Linear predictor for a single MCMC sample
mcmc_sample_lp <- function(sample, object)
{
    D <- as.vector(dist(matrix(sample[-object$beta_idx], ncol=object$d)))
    object$X %*% matrix(coef(object), ncol=1) - D
}

predict.lsmfit <- function(object, type="link", ...)
{
    if (object$method == "MLE") {
        D <- as.vector(dist(locations(object)))
        pred <- object$X %*% matrix(coef(object), ncol=1) - D
    } else if (object$method == "MH") {
        pred <- apply(object$estimate$samples, 1,
                      function(s) mcmc_sample_lp(s, object))
    } else {
        stop("unknown estimation method")
    }

    if (type == "response") {
        if (object$family == "bernoulli") {
            pred <- plogis(pred)
        } else {
            stop("unsupported family")
        }
    }

    pred
}
