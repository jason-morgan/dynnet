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
        est <- summary(model$estimate$samples)$statistics[idx+1,"Mean"]
    }

    est
}

summary.lsmfit <- function(object, ...)
{
    idx <- object$beta_idx
    est <- coef(object)
    se  <- sqrt(diag(vcov(object))[idx])

    tbl <- data.frame("Estimate"=est, "SE"=se)
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
