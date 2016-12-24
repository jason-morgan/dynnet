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

    names(est) <- "(Intercept)"
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

    rownames(tbl) <- "(Intercept)"
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
