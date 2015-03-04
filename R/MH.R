lsm <- function(network, start_b, start_Z, ref_idx, family, iter=10000)
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

    smry    <- list(b_accepted=0, Z_accepted=0)
    current <- list(b=start_b, Z=start_Z, posterior=samples[1,p_idx], smry=smry)

    for (i in 2:iter) {
        if (i %% 1000 == 0) cat("Iter:", i, "\n")

        current <- MH_iter(current$b,
                           current$Z,
                           ref_idx,
                           y, X,
                           current$posterior,
                           current$smry,
                           log_posterior)

        samples[i,b_idx] <- current$b
        samples[i,Z_idx] <- as.vector(current$Z)
        samples[i,p_idx] <- current$posterior
    }

    cat("\nb accepted:", current$smry$b_accepted / iter,
        "\nZ accepted:", current$smry$Z_accepted / iter, "\n")

    ## mcmc(data=samples, start=0, end=iter, thin=thin)
    list(samples, idx=list(ref=ref_idx, b=b_idx, Z=Z_idx, lp=p_idx))
}

MH_iter <- function(b, Z, ref_idx, y, X, cur_post, smry, log_posterior)
{
    b_prop  <- list(b=propose_b(b), Z=Z)
    current <- list(b=b, Z=Z, posterior=cur_post)
    b_step  <- MH_step(b_prop, current, ref_idx, y, X, log_posterior)

    Z_prop  <- list(b=b, Z=propose_Z(Z, ref_idx, i=5))
    current <- list(b=b_step$b, Z=Z, posterior=b_step$posterior)
    Z_step  <- MH_step(Z_prop, current, ref_idx, y, X, log_posterior)

    ## Z_prop  <- list(b=propose_b(b), Z=propose_Z(Z, ref_idx))
    ## current <- list(b=b, Z=Z, posterior=cur_post)
    ## Z_step  <- MH_step(Z_prop, current, ref_idx, y, X, log_posterior)

    smry$b_accepted <- smry$b_accepted + b_step$accept
    smry$Z_accepted <- smry$Z_accepted + Z_step$accept

    list(b=Z_step$b, Z=Z_step$Z, posterior=Z_step$posterior, smry=smry)
}

MH_step <- function(proposal, current, ref_idx, y, X, log_posterior)
{
    new_post <- log_posterior(y, proposal$b, proposal$Z, X, ref_idx)
    p <- exp(new_post - current$posterior)

    if (runif(1) < p) {
        list(b=proposal$b, Z=proposal$Z, posterior=new_post, accept=1)
    } else {
        list(b=current$b, Z=current$Z, posterior=current$posterior, accept=0)
    }
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

propose_b <- function(b)
{
    b + rnorm(length(b), mean=0, sd=0.5)
}

propose_Z <- function(Z, ref_idx, i=nrow(Z)-length(ref_idx))
{
    idx     <- sample(seq_len(nrow(Z))[-ref_idx], i)
    ## S       <- (length(Z) - (ncol(Z) * length(ref_idx))) - 2
    S       <- i * ncol(Z)

    ## idx     <- seq_len(nrow(Z))[-ref_idx]
    ## S       <- length(idx) * ncol(Z)

    e       <- mvtnorm::rmvnorm(i, mean=rep(0, ncol(Z)),
                                sigma=diag(1/S, nrow=ncol(Z)))
    Z[idx,] <- Z[idx,] + e
    Z
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

log_prior <- function(b, Z, ref_idx, b_sd=5, Z_sd=5)
{
    b_sigma <- diag(b_sd, nrow=length(b))
    Z_sigma <- diag(Z_sd, nrow=ncol(Z))

    B_prior <- sum(mvtnorm::dmvnorm(b, sigma=b_sigma, log=TRUE))
    Z_prior <- sum(mvtnorm::dmvnorm(Z[-ref_idx,], sigma=Z_sigma, log=TRUE))

    B_prior + Z_prior
}
