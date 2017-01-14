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
    ## shift idx +1 to account for log probability
    idx <- object$beta_idx + 1
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
        tbl <- data.frame("Estimate"=est, "Std. Error"=se)
    } else if (object$method == "MH") {
        if (length(idx) == 1)
            se <- sd(object$estimate$samples[, idx])
        else
            se <- apply(object$estimate$samples[, idx], 2, sd)

        tbl <- data.frame("Posterior Mean"=est, "Std. Dev."=se)
    }

    rownames(tbl) <- object$beta_names

    cat("Latent Space Model\n")
    cat("Call:\n")
    print(object$call)

    cat("Estimation method: ", object$method, "\n")

    cat("\nCoefficients on exogenous predictors:\n")
    print(tbl)
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
##' model. For a model fit via MLE, this is the MLE estimates. For models with
##' with the Metropolis-Hastings algorithm, the posterior means for the
##' locations (possibly transformed via a Procrustes step, if that was the
##' method of identification specified). Note: the posterior means should not be
##' used for inference as they can distort the appearance of vertices in the the
##' latent space. Use the posterior samples instead.
##' @title Estimated Locations from Latent Space Model
##' @param object A \code{\link{lsmfit}} object.
##' @param ... Addition parameters to pass.
##' @return A matrix of estimated locations. See Details.
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
        est <- cbind(est, mean(est))

    rownames(est) <- vertex_attr(object$graph, "name")
    est
}

## Linear predictor for a single MCMC sample, includes reference units when the
## model had them.
mcmc_sample_lp <- function(sample, model)
{
    theta <- decompose_theta(sample, model$beta_idx, model$d)

    if (!is.null(model$ref)) {
        Z <- insert_ref(T$Z, model$ref, model$d)
    } else {
        Z <- T$Z
    }

    D <- as.vector(dist(Z))
    model$X %*% matrix(theta$beta, ncol=1) - D
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

##' Print dynnet network object
##'
##' Print dynnet network object
##' @title Print dynnet Network Object
##' @param network dynnet network object.
##' @return Prints.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
print <- function(object)
{
    UseMethod("print")
}

##' @export
##' @rdname print
print.dynnet <- function(object)
{
    cat("DYNNET\n")
    cat("  |- Periods:", periods(object), "\n")
    cat("  |- Nodes:  ", do.call(c, vcount(object)), "\n")
}

##' @export
##' @rdname print
print.lsmfit <- function(object)
{
    cat("Latent Space Model\n")
    cat("  |- Call:")
    print(object$call)

    cat("  |- Estimation method: ", object$method, "\n")
    cat("  |- Distance metric: ", object$dist_metric, "\n")
}
