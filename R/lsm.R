lsm <- function(network, k=1, period=1, ref=NULL, family="bernoulli",
                start=start_random, method="MLE", seed=NULL, verbose=TRUE,
                control=list(MCMC.samplesize=2^10, MCMC.interval=10))
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

    ## Starting values
    pos <- start(graph, ref, k)
    beta0 <- 1
    theta <- c(beta0, as.vector(pos[-ref$idx,]))

    model <- structure(list(edges=edges, period=1, k=k, ref=ref,
                            family=family,
                            dropped=which(deg == 0), seed=seed),
                       class="dynnetlsm")

    if (method == "MLE")
        est <- lsm_MLE(theta, model, control=control)
    else if (method == "MH")
        est <- lsm_MH(theta, model, control=control)
    else
        est <- NULL

    model$graph <- graph
    model$estimate <- est

    model
}

lsm_MLE <- function(theta, model, control)
{
    est <- optim(theta, calc_likelihood, model=model, method="L-BFGS-B",
                 control=list(trace=1, fnscale=-1, maxit=500))

    est
}

insert_ref <- function(pos, ref, k)
{
    n <- nrow(pos) + length(ref$idx)

    .C_insert_ref(ref$idx, ref$pos, (1:n)[-ref$idx], pos)
}

calc_likelihood <- function(theta, model=NULL)
{
    beta0 <- theta[1]
    pos <- insert_ref(matrix(theta[-1], ncol=model$k), model$ref, model$k)
    llik_logit(model$edges, beta0 - .C_dist_euclidean(pos))
}

start_random <- function(graph, ref, k)
{
    n <- igraph::vcount(graph)
    pos <- matrix(runif(n*k, min=-0.5, max=0.5), ncol=k, nrow=n)
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
