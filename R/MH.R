lsm_MH <- function(network, start_b, start_Z, ref_idx, family, iter=10000)
{
    y <- as.matrix(network$adj[[1]])
    y <- y[lower.tri(y)]
    X <- network$X[[1]]

    log_posterior <- mk_log_posterior(family)

    k       <- ncol(start_Z)
    b_len   <- length(start_b)
    Z_len   <- length(start_Z)

    samples <- matrix(0, ncol=(1 + b_len + Z_len), nrow=iter)
    colnames(samples) <- c(paste0("b", 0:(b_len-1)),
                           paste0("Z", 1:Z_len),
                           "lp")

    b_idx <- 1:b_len
    Z_idx <- (b_len+1):(b_len+Z_len)
    p_idx <- ncol(samples)

    samples[1,b_idx] <- start_b
    samples[1,Z_idx] <- as.vector(start_Z)
    samples[1,p_idx] <- log_posterior(y, start_b, start_Z, X, ref_idx)

    model <- list(y=y, X=X, ref_idx=ref_idx, family=family,
                  log_posterior=log_posterior)
    state <- list(theta=list(b=start_b, Z=start_Z, posterior=samples[1,p_idx]),
                  smry=list(b_accepted=0, b_iter=0, Z_accepted=0, Z_iter=0))

    for (i in 2:iter) {
        if (i %% 1000 == 0) cat("Iter:", i, "\n")

        state <- MCMC_step_Z(state, model)
        state <- MCMC_step_b(state, model)

        samples[i,b_idx] <- state$theta$b
        samples[i,Z_idx] <- as.vector(state$theta$Z)
        samples[i,p_idx] <- state$theta$posterior
    }

    cat("\nb accepted:", state$smry$b_accepted / state$smry$b_iter,
        "\nZ accepted:", state$smry$Z_accepted / state$smry$Z_iter, "\n")

    list(samples, idx=list(ref=ref_idx, b=b_idx, Z=Z_idx, lp=p_idx))
}

MCMC_step_Z <- function(state, model)
{
    ## idx <- seq_len(nrow(state$theta$Z))[-model$ref_idx]
    ## order_idx <- sample(idx, length(idx), replace=TRUE)

    ## Update the latent locations of each node separately.
    ## for (i in order_idx) {
    ##     proposal <- propose_Z(i, state)
    ##     update   <- MH_update(proposal, state, model)
    ##     state    <- update$state
    ##     state$smry$Z_accepted <- state$smry$Z_accepted + update$accept
    ##     state$smry$Z_iter <- state$smry$Z_iter + 1
    ## }

    proposal <- propose_Z(model$ref_idx, state)
    update   <- MH_update(proposal, state, model)
    state    <- update$state
    state$smry$Z_accepted <- state$smry$Z_accepted + update$accept
    state$smry$Z_iter <- state$smry$Z_iter + 1

    state
}

MCMC_step_b <- function(state, model)
{
    proposal <- propose_b(state)
    update   <- MH_update(proposal, state, model)
    state    <- update$state
    state$smry$b_accepted <- state$smry$b_accepted + update$accept
    state$smry$b_iter <- state$smry$b_iter + 1

    state
}

MH_update <- function(proposal, state, model)
{
    new_post <- model$log_posterior(model$y, proposal$theta$b, proposal$theta$Z,
                                    model$X, model$ref_idx)

    p <- exp(new_post - state$theta$posterior)
    accept <- 0

    if (runif(1) < p) {
        state$theta <- proposal$theta
        state$theta$posterior <- new_post
        accept <- 1
    }

    list(state=state, accept=accept)
}

mk_log_posterior <- function(family="logit")
{
    llik <- llik_fn(family)

    function(y, b, Z, X, ref_idx)
    {
        d  <- as_distance_vector(Z)
        lp <- (X %*% b) - d
        llik(y, lp) + log_prior(b, Z, ref_idx)
    }
}

propose_b <- function(state)
{
    scale <- 1 / (length(state$theta$b) + 1)
    e <- rnorm(length(state$theta$b), mean=0, sd=scale)
    state$theta$b <- state$theta$b + e
    state
}

propose_Z <- function(ref_idx, state)
{
    n <- nrow(state$theta$Z[-ref_idx,])
    N <- n * ncol(state$theta$Z)
    scale <- 1 / N

    e <- mvtnorm::rmvnorm(n,
                          mean=rep(0, ncol(state$theta$Z)),
                          sigma=diag(ncol(state$theta$Z))*scale)

    state$theta$Z[-ref_idx,] <- state$theta$Z[-ref_idx,] + e
    state
}

propose_one_Z <- function(idx, state)
{
    scale <- 1 / (ncol(state$theta$Z) - 1)
    e <- rnorm(ncol(state$theta$Z), mean=0, sd=scale)

    state$theta$Z[idx,] <- state$theta$Z[idx,] + e
    state
}

llik_logit <- function(y, lp)
{
    sum(dbinom(y, 1, plogis(lp), log=TRUE))
}

llik_poisson <- function(y, lp)
{
    sum(dpois(y, 1, exp(lp), log=TRUE))
}

llik_fn <- function(family)
{
    switch(family,
           "logit"   = llik_logit,
           "poisson" = llik_poisson,
           stop("unknown family"))
}

log_prior <- function(b, Z, ref_idx, b_sd=10, Z_sd=5)
{
    b_sigma <- diag(b_sd, nrow=length(b))
    Z_sigma <- diag(Z_sd, nrow=ncol(Z))

    B_prior <- sum(mvtnorm::dmvnorm(b, sigma=b_sigma, log=TRUE))
    Z_prior <- sum(mvtnorm::dmvnorm(Z[-ref_idx,], sigma=Z_sigma, log=TRUE))

    B_prior + Z_prior
}
