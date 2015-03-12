## lsm_MLE <- function(network, k=1, start=NULL, family="logit")
lsm_MLE <- function(model)
{
    llik <- mk_log_likelihood(model$family, k=model$k)

    ## Hacky
    y <- get_adjacency(model$network, period=model$period)
    y <- as.matrix(Matrix::tril(y))
    y <- y[lower.tri(y)]

    X     <- model$X
    b_idx <- 1:ncol(X)
    Z_idx <- (length(b_idx)+1):(length(model$start))
    template <- matrix(model$start[Z_idx], ncol=model$k)

    theta <- c(model$start[b_idx],
               as.vector(matrix(model$start[Z_idx], ncol=model$k)[-ref$idx,]))

    Z_idx <- (length(b_idx)+1):(length(theta))

    tol <- sqrt(.Machine$double.eps)

    ## Initial estimates of Z
    cur <- optim(theta[Z_idx], MLE_Z_update, llik=llik, theta=theta, y=y, X=X,
                 b_idx=b_idx, Z_idx=Z_idx, ref_idx=ref$idx, template=template,
                 ## method="BFGS",
                 control=list(trace=0, fnscale=-1))

    theta[Z_idx] <- cur$par
    ## theta <- center_Z(theta, b_idx, Z_idx, k)

    old_lik <- 0
    new_lik <- 1
    iter <- 0
    MAXITER <- 100
    while (MAXITER > iter && (abs(new_lik - old_lik) > tol)) {
        old_lik <- new_lik
        cur <- optim(theta[b_idx], MLE_b_update, llik=llik, theta=theta, y=y,
                     X=X, b_idx=b_idx, Z_idx=Z_idx, ref_idx=ref$idx, template=template,
                     method="Brent", lower=-100, upper=100,
                     control=list(trace=0, fnscale=-1))
        theta[b_idx] <- cur$par

        cur <- optim(theta[Z_idx], MLE_Z_update, llik=llik, theta=theta, y=y, X=X,
                     b_idx=b_idx, Z_idx=Z_idx, ref_idx=ref$idx, template=template,
                     control=list(trace=0, fnscale=-1))
        theta[Z_idx] <- cur$par
        ## theta <- center_Z(theta, b_idx, Z_idx, k)

        if (iter %% 2 == 0) {
            cat("Iteration:", iter, "\n")
            cat("  beta:", theta[b_idx], "\n")
            cat("  log-likelihood:", cur$value, "\n")
        }
        iter <- iter + 1
        new_lik <- cur$value
    }

    model$b_idx <- b_idx
    model$Z_idx <- Z_idx
    list(theta=theta, last=cur, model=model)
}

center_Z <- function(theta, b_idx, Z_idx, k)
{
    Z_est <- scale(matrix(theta[Z_idx], ncol=k), scale=FALSE)
    c(theta[b_idx], as.vector(Z_est))
}

MLE_Z_update <- function(Z, llik=NULL, theta=NULL, y=NULL, X=NULL, b_idx=NULL,
                         Z_idx=NULL, ref_idx=NULL, template=NULL)
{
    theta[Z_idx] <- Z
    llik(theta, y=y, X=X, b_idx=b_idx, Z_idx=Z_idx, ref_idx=ref_idx,
         template=template)
}

MLE_b_update <- function(b, llik=NULL, theta=NULL, y=NULL, X=NULL, b_idx=NULL,
                         Z_idx=NULL, ref_idx=NULL, template=NULL)
{
    theta[b_idx] <- b
    llik(theta, y=y, X=X, b_idx=b_idx, Z_idx=Z_idx, ref_idx=ref_idx,
         template=template)
}

make_Z <- function(theta, Z_idx, template, ref_idx, k)
{
    template[-ref_idx,] <- matrix(theta[Z_idx], ncol=k)
    template
}

mk_log_likelihood <- function(family, k)
{
    llik <- llik_fn(family)

    function(theta, y=NULL, X=NULL, b_idx=NULL, Z_idx=NULL, ref_idx=NULL,
             template=NULL)
    {
        b  <- theta[b_idx]
        Z  <- make_Z(theta, Z_idx, template, ref_idx, k)
        d  <- as_distance_vector(Z)
        lp <- (X %*% b) - d
        llik(y, lp)
    }
}
