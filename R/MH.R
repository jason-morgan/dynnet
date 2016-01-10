lsm_MH <- function(theta, model, control=list(MCMC.samplesize=2^10,
                                              MCMC.interval=10))
{
    beta0 <- theta[1]
    Z     <- matrix(theta[-1], ncol=model$k)
    pos   <- insert_ref(Z, model$ref, model$k)
    D     <- as.vector(.C_dist_euclidean(pos))
    posterior <- .C_log_posterior_logit(model$edges, beta0 - D, beta0, Z)

    state <- list(edges=model$edges,
                  beta0=beta0,
                  Z=Z,
                  ref=model$ref,
                  k=model$k,
                  posterior=posterior,
                  iter=0,
                  accept=0)

    BURN <- control$MCMC.samplesize * 100
    ## BURN <- 2^23

    ## burnin
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
    record <- matrix(0, ncol=(length(beta0) + length(Z) + 1),
                     nrow=control$MCMC.samplesize)
    colnames(record) <- c("posterior", "beta0", paste0("z", 1:length(Z)))

    i <- 1
    cat("SAMPLING (", control$MCMC.samplesize, "): ", sep="")
    while (i <= control$MCMC.samplesize) {
        state <- MH_sampler(state, control$MCMC.interval)
        record[i,] <- c(state$posterior, state$beta0, as.vector(state$Z))
        if (i %% 100 == 0) cat(i, " ")
        i <-  i + 1
    }
    cat("done\n")

    cat("SAMPLING Acceptance rate:", state$accept / state$iter, "\n")

    list(state=state, record=coda::mcmc(record))
}

MH_sampler <- function(state, n_iter)
{
    for (i in 1:n_iter) {
        proposal <- list(Z     = .C_propose_Z(state$Z),
                         beta0 = .C_propose_beta0(state$beta0))
        state    <- MH_update(proposal, state)
    }

    state
}

MH_update <- function(proposal, state)
{
    state$iter <- state$iter + 1

    pos <- insert_ref(proposal$Z, state$ref, state$k)
    D   <- .C_dist_euclidean(pos)
    proposal$posterior <- .C_log_posterior_logit(state$edges,
                                                 proposal$beta0 - D,
                                                 proposal$beta0,
                                                 proposal$Z)

    p <- exp(proposal$posterior - state$posterior)

    if (runif(1) < p) {
        state$beta0  <- proposal$beta0
        state$Z      <- proposal$Z
        state$posterior <- proposal$posterior
        state$accept <- state$accept + 1
    }

    state
}
