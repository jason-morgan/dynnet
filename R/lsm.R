lsm <- function(network, k=1, period=1, ref=NULL,
                family="bernoulli", method="MLE", fit=TRUE, seed=NULL,
                verbose=TRUE)
{
    if (!is.null(seed))
        set.seed(seed)

    if (is.null(ref))
        stop("lsm: reference node information not provided")

    ## Remove isolates; their distance from all other nodes is infinite
    graph <- get_graph(network, period=period)
    deg   <- igraph::degree(graph)
    graph <- igraph::delete.vertices(graph, deg == 0)

    if (sum(deg == 0) > 0)
        warning("lsm: isolate(s) removed from network", call.=FALSE)

    ## This will need to be modified for directed networks
    edges <- extract_edges(graph)

    ## Set up latent positions matrix
    pos <- starting_positions(graph, ref, k)

    beta0 <- qlogis(igraph::edge_density(graph))
    print(beta0)
    theta <- c(beta0, as.vector(pos[-ref$idx,]))

    model <- structure(list(edges=edges, period=1, k=k, ref=ref,
                            family=family,
                            dropped=which(deg == 0), seed=seed),
                       class="dynnetlsm")

    est <- optim(theta, calc_likelihood, model=model, method="L-BFGS-B",
                 control=list(trace=1, fnscale=-1, maxit=250))

    ## est <- optim(est$par, calc_likelihood, model=model, method="SANN",
    ##              control=list(trace=1, fnscale=-1, maxit=500000, tmax=50))

    model$graph <- graph
    model$estimate <- est
    model
}

insert_ref <- function(pos, ref, k)
{
    n <- nrow(pos) + length(ref$idx)
    new <- matrix(rnorm(n*k), ncol=k, nrow=n)
    new[(1:n)[-ref$idx],] <- pos
    new[ref$idx,] <- ref$pos
    new
}

calc_likelihood <- function(theta, model=NULL)
{
    beta0 <- theta[1]
    pos <- insert_ref(matrix(theta[-1], ncol=model$k), model$ref, model$k)
    llik_logit(model$edges, beta0 - as.matrix(dist(pos)))
}

starting_positions <- function(graph, ref, k)
{
    n <- igraph::vcount(graph)
    ## For the moment, set random positions
    ## pos <- matrix(rnorm(n*k), ncol=k, nrow=n)
    pos <- matrix(runif(n*k, min=-0.5, max=0.5), ncol=k, nrow=n)
    ## pos <- matrix(0, ncol=k, nrow=n)
    pos[ref$idx,] <- ref$pos
    pos
}

extract_edges <- function(graph)
{
    y <- as.matrix(as_adj(graph))
    y[lower.tri(y)]
}

llik_fn <- function(family)
{
    switch(family,
           "bernoulli" = llik_logit,
           "poisson"   = llik_poisson,
           stop("unknown family"))
}

plot.dynnetlsm <- function(model)
{
    plot(rbind(model$ref$pos, model$estimate$par))
}

log_posterior_logit <- function(llik_fn, y, lp)
{
    llik(y, lp) + log_prior_b(b) + log_prior_Z(Z, ref_idx)
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
