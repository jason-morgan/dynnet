##' Static latent space model for static networks.
##'
##' Estimate a latent space model for static networks as introduced by HRH
##' (2002).
##' @title Static Latent Space Model
##' @param network \code{dynnet} network object.
##' @param k Integer. Number of dimensions for the latent space.
##' @param period Integer. Period to select from the network object.
##' @param ref List. Reference vertex information.
##' @param family String specifying the model family. Only \code{bernoulli} is
##'     currently supported.
##' @param start Function to use to get the starting values for the latent
##'     locations of each vertex.
##' @param method String specifying the estimation method. Either \code{MLE}
##'     (maximum likelihood) or \code{MH} (a Metropolis Hastings algorithm).
##' @param seed Random seed.
##' @param verbose Boolean.
##' @param control List of control parameters to be used for
##' @return \code{lsmfit} model object.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
lsm <- function(network, k=1, period=1, ref=NULL, family="bernoulli",
                start=start_random, method="MLE", seed=NULL, verbose=TRUE,
                control=control.lsm())
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
    alpha <- 1
    theta <- c(alpha, as.vector(pos[-ref$idx,]))

    model <- structure(list(edges=edges, period=1, k=k, ref=ref,
                            family=family, beta_idx=1,
                            dropped=which(deg == 0), seed=seed,
                            verbose=verbose),
                       class="lsmfit")

    if (method == "MLE") {
        est <- lsm_MLE(theta, model, control=control)
    } else if (method == "MH") {
        theta <- lsm_MLE(theta, model, control=control)$par
        est   <- lsm_MH(theta, model, control=control)

        ## remove fixed values
        rm_idx <- which(apply(est$samples, 2, var) == 0)
        est$samples <- coda::mcmc(est$samples[,-rm_idx])
    } else {
        est <- NULL
    }

    model$method <- method
    model$graph <- graph
    model$estimate <- est

    model
}

lsm_MLE <- function(theta, model, control)
{
    est <- optim(theta, calc_likelihood, model=model, method="L-BFGS-B",
                 control=list(trace=model$verbose, fnscale=-1, maxit=500),
                 hessian=TRUE)

    est
}

insert_ref <- function(pos, ref, k)
{
    n <- nrow(pos) + length(ref$idx)

    .C_insert_ref(ref$idx, ref$pos, (1:n)[-ref$idx], pos)
}

calc_likelihood <- function(theta, model=NULL)
{
    alpha <- theta[model$beta_idx]
    pos <- insert_ref(matrix(theta[-model$beta_idx], ncol=model$k),
                      model$ref, model$k)
    .C_llik_logit(model$edges, alpha - .C_dist_euclidean(pos))
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
