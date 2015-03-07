lsm_MH <- function(network, start_b, start_Z, ref_idx, family,
                   burnin=10000, nsamples=10000)
{
    y <- as.matrix(network$adj[[1]])
    y <- y[lower.tri(y)]
    X <- network$X[[1]]

    log_posterior <- mk_log_posterior(family)

    k       <- ncol(start_Z)
    b_len   <- length(start_b)
    Z_len   <- length(start_Z)

    samples <- matrix(0, ncol=(1 + b_len + Z_len), nrow=nsamples)
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

    ## Burnin stage
    cat("Beginning burnin stage: ")
    for (i in 2:burnin) {
        if (i %% 1000 == 0) cat(i, "...", sep="")

        state <- MCMC_step_Z(state, model)
        state <- MCMC_step_b(state, model)
    }

    cat("\nBurnin stage acceptance rates:")
    cat("\n  b accepted:", state$smry$b_accepted / state$smry$b_iter,
        "\n  Z accepted:", state$smry$Z_accepted / state$smry$Z_iter, "\n")


    ## Sampling stage
    state$smry$b_accepted <- 0
    state$smry$b_iter     <- 0
    state$smry$Z_accepted <- 0
    state$smry$Z_iter     <- 0

    cat("\nBeginning sampling stage: ")
    for (i in 1:nsamples) {
        if (i %% 1000 == 0) cat(i, "...", sep="")

        state <- MCMC_step_Z(state, model)
        state <- MCMC_step_b(state, model)

        samples[i,b_idx] <- state$theta$b
        samples[i,Z_idx] <- as.vector(state$theta$Z)
        samples[i,p_idx] <- state$theta$posterior
    }

    cat("\nSampling stage acceptance rates:")
    cat("\n  b accepted:", state$smry$b_accepted / state$smry$b_iter,
        "\n  Z accepted:", state$smry$Z_accepted / state$smry$Z_iter, "\n")

    list(samples, idx=list(ref=ref_idx, b=b_idx, Z=Z_idx, lp=p_idx))
}

MCMC_step_Z <- function(state, model)
{
    idx <- seq_len(nrow(state$theta$Z))[-model$ref_idx]
    order_idx <- sample(idx, length(idx), replace=TRUE)

    ## Update the latent locations of each node separately.
    for (i in order_idx) {
        proposal <- propose_one_Z(i, state)
        update   <- MH_update(proposal, state, model)
        state    <- update$state
        state$smry$Z_accepted <- state$smry$Z_accepted + update$accept
        state$smry$Z_iter <- state$smry$Z_iter + 1
    }

    ## proposal <- propose_Z(model$ref_idx, state)
    ## update   <- MH_update(proposal, state, model)
    ## state    <- update$state
    ## state$smry$Z_accepted <- state$smry$Z_accepted + update$accept
    ## state$smry$Z_iter <- state$smry$Z_iter + 1

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
        llik(y, lp) + log_prior_b(b) + log_prior_Z(Z, ref_idx)
    }
}

propose_b <- function(state)
{
    ## scale <- 1 / (length(state$theta$b) + 1)
    ## e <- rnorm(length(state$theta$b), mean=0, sd=scale)
    e <- runif(1, min=-2, max=2)
    state$theta$b <- abs(state$theta$b + e)
    state
}

propose_Z <- function(ref_idx, state)
{
    n <- nrow(state$theta$Z[-ref_idx,])
    N <- n * ncol(state$theta$Z)
    scale <- 1 / sqrt(N)

    e <- mvtnorm::rmvnorm(n,
                          mean=rep(0, ncol(state$theta$Z)),
                          sigma=diag(ncol(state$theta$Z))*scale)

    state$theta$Z[-ref_idx,] <- state$theta$Z[-ref_idx,] + e
    state
}

propose_one_Z <- function(idx, state)
{
    ## scale <- 1 / (ncol(state$theta$Z) - 1)
    scale <- 1.5
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

log_prior_b <- function(b, b_sd=100)
{
    b_sigma <- diag(b_sd, nrow=length(b))
    B_prior <- sum(mvtnorm::dmvnorm(b, sigma=b_sigma, log=TRUE))

    ## Matches the prior used by HRH (2002)
    ## B_prior <- dgamma(b, 1, scale=1, log=TRUE)
    B_prior
}

log_prior_Z <- function(Z, ref_idx, Z_sd=100)
{
    Z_sigma <- diag(Z_sd, nrow=ncol(Z))
    Z_prior <- sum(mvtnorm::dmvnorm(Z[-ref_idx,], sigma=Z_sigma, log=TRUE))

    Z_prior
}
