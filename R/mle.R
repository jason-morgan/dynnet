lsm_MLE <- function(network, k=1, start=NULL, family="logit")
{
    llik <- mk_log_likelihood(family, k=k)

    y <- network$adj[[1]][lower.tri(network$adj[[1]])]
    X <- network$X[[1]]
    b_idx <- 1:ncol(X)
    Z_idx <- (length(b_idx)+1):(nrow(network$adj[[1]])*k+length(b_idx))

    if (is.null(start))
        theta <- c(rep(0, length(b_idx)), rnorm(length(Z_idx)))
    else
        theta <- center_Z(start, b_idx, Z_idx, k)

    tol <- sqrt(.Machine$double.eps)

    ## Initial estimates of Z
    cur <- optim(theta[Z_idx], MLE_Z_update, llik=llik, theta=theta, y=y, X=X,
                 b_idx=b_idx, Z_idx=Z_idx, method="BFGS",
                 control=list(trace=0, fnscale=-1))
    theta[Z_idx] <- cur$par
    theta <- center_Z(theta, b_idx, Z_idx, k)

    cur_lik <- 0
    iter <- 0
    MAXITER <- 100
    while (iter < MAXITER || (abs(cur$value - cur_lik) > tol )) {
        cur <- optim(theta[b_idx], MLE_b_update, llik=llik, theta=theta, y=y,
                     X=X, b_idx=b_idx, Z_idx=Z_idx, method="Brent", lower=-100,
                     upper=100, control=list(trace=0, fnscale=-1))
        theta[b_idx] <- cur$par

        cur <- optim(theta[Z_idx], MLE_Z_update, llik=llik, theta=theta, y=y, X=X,
                     b_idx=b_idx, Z_idx=Z_idx,
                     control=list(trace=0, fnscale=-1))
        theta[Z_idx] <- cur$par
        theta <- center_Z(theta, b_idx, Z_idx, k)

        cat("Log-likelihood:", cur$value, "\n")
        iter <- iter + 1
        cur_lik <- cur$value
    }

    list(theta=theta, last=cur)
}

center_Z <- function(theta, b_idx, Z_idx, k)
{
    Z_est <- scale(matrix(theta[Z_idx], ncol=k), scale=FALSE)
    c(theta[b_idx], as.vector(Z_est))
}

MLE_Z_update <- function(Z, llik=NULL, theta=NULL, y=NULL, X=NULL, b_idx=NULL,
                         Z_idx=NULL)
{
    theta[Z_idx] <- Z
    llik(theta, y=y, X=X, b_idx=b_idx, Z_idx=Z_idx)
}

MLE_b_update <- function(b, llik=NULL, theta=NULL, y=NULL, X=NULL, b_idx=NULL,
                         Z_idx=NULL)
{
    theta[b_idx] <- b
    llik(theta, y=y, X=X, b_idx=b_idx, Z_idx=Z_idx)
}

mk_log_likelihood <- function(family, k)
{
    llik <- llik_fn(family)

    function(theta, y=NULL, X=NULL, b_idx=NULL, Z_idx=NULL)
    {
        b  <- theta[b_idx]
        Z  <- matrix(theta[Z_idx], ncol=k)
        d  <- as_distance_vector(Z)
        lp <- (X %*% b) - d
        llik(y, lp)
    }
}
