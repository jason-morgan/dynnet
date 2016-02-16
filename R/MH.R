lsm_MH <- function(theta, model, control=list(MCMC.samplesize=2^10,
                                              MCMC.interval=10))
{
    alpha <- theta[model$beta_idx]
    Z     <- matrix(theta[-model$beta_idx], ncol=model$k)
    pos   <- insert_ref(Z, model$ref, model$k)
    D     <- as.vector(.C_dist_euclidean(pos))
    posterior <- .C_log_posterior_logit(model$edges, alpha - D, alpha, Z)

    state <- list(edges=model$edges,
                  alpha=alpha,
                  Z=Z,
                  ref=model$ref,
                  k=model$k,
                  posterior=posterior,
                  iter=0,
                  accept=0)

    ## burnin
    BURN <- control$MCMC.burnin

    i <- 1
    cat("BURNIN (", BURN, "): ", sep="")
    while (i <= 10) {
        state <- MH_sampler(state, BURN/10)
        i <- i + 1
        cat(i * BURN/10, " ")
    }
    cat("done\n")

    cat("BURNIN Acceptance rate:", state$accept / state$iter, "\n")
    state$iter <- 0
    state$accept <- 0

    ## sampling
    samples <- matrix(0, ncol=(length(alpha) + length(Z) + 1),
                     nrow=control$MCMC.samplesize)
    colnames(samples) <- c("posterior", "alpha", paste0("z", 1:length(Z)))

    i <- 1
    cat("SAMPLING (", control$MCMC.samplesize, "): ", sep="")
    while (i <= control$MCMC.samplesize) {
        state <- MH_sampler(state, control$MCMC.interval)
        samples[i,] <- c(state$posterior, state$alpha, as.vector(state$Z))
        if (i %% 100 == 0) cat(i, " ")
        i <-  i + 1
    }
    cat("done\n")

    cat("SAMPLING Acceptance rate:", state$accept / state$iter, "\n")

    list(state=state, samples=coda::mcmc(samples))
}

MH_sampler <- function(state, n_iter)
{
    for (i in 1:n_iter) {
        proposal <- list(Z     = .C_propose_Z(state$Z),
                         alpha = .C_propose_alpha(state$alpha))
        state    <- MH_update(proposal, state)
    }

    state
}

MH_update <- function(proposal, state)
{
    state$iter <- state$iter + 1

    n <- nrow(proposal$Z) + length(state$ref$idx)
    est_idx <- (1:n)[-state$ref$idx]

    pos <- .C_insert_ref(state$ref$idx, state$ref$pos, est_idx, proposal$Z)
    D   <- .C_dist_euclidean(pos)
    proposal$posterior <- .C_log_posterior_logit(state$edges,
                                                 proposal$alpha - D,
                                                 proposal$alpha,
                                                 proposal$Z)

    p <- exp(proposal$posterior - state$posterior)

    if (runif(1) < p) {
        state$alpha  <- proposal$alpha
        state$Z      <- proposal$Z
        state$posterior <- proposal$posterior
        state$accept <- state$accept + 1
    }

    state
}

lsm_MH_C <- function(theta, model, control=list(MCMC.samplesize=2^10,
                                                MCMC.interval=10))
{
    alpha <- theta[1]
    Z     <- insert_ref(matrix(theta[-1], ncol=model$k), model$ref, model$k)
    Z_idx <- (1:nrow(Z))[-model$ref$idx]
    X <- matrix(0.0)
    beta <- c(0)

    .C_lsm_MH(model$edges, X, Z_idx, model$k,
              control$MCMC.burnin, control$MCMC.samplesize,
              control$MCMC.interval, alpha, beta, Z)
}
