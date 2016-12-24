##' Static latent space model for static networks.
##'
##' Estimate a latent space model for static networks as introduced by HRH
##' (2002).
##' @title Static Latent Space Model
##' @param formula A \code{\link{formula}} object of the form \code{network ~
##'     terms}, where \code{code} is a \code{dynnet} network object and
##'     \code{terms} are network terms to be included in the model.
##' @param d Integer. Number of dimensions for the latent space.
##' @param period Integer. Period to select from the network object.
##' @param ref List. Reference vertex information. If \code{NULL}, then
##'     reference nodes will not be used for identification.
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
lsm <- function(formula, d=1, period=1, ref=NULL, family="bernoulli",
                start=start_random, method="MLE", seed=NULL, verbose=TRUE,
                control=control.lsm())
{
    if (!is.null(seed))
        set.seed(seed)

    fc <- match.call(expand.dots=TRUE)
    ft <- terms(formula, specials=c("absdiff"))
    network <- get(all.vars(ft)[[attr(ft, "response")]])
    graph <- get_graph(network, period=period)

    ## Remove isolates; their distance from all other nodes is infinite
    deg   <- igraph::degree(graph)
    graph <- igraph::delete.vertices(graph, deg == 0)

    if (sum(deg == 0) > 0)
        warning("lsm: isolate(s) removed from network", call.=FALSE,
                immediate.=TRUE)

    ## Build model matrix and extract edge values
    dta <- build_model_data(ft, graph)
    edges <- dta$edges
    X <- dta$x

    ## Starting values
    pos  <- start(graph, ref, d)
    beta <- rep(0, ncol(X))
    if (!is.null(ref))
        theta <- c(beta, as.vector(pos[-ref$idx,]))
    else
        theta <- c(beta, as.vector(pos))

    ## Set node names
    nodes <- igraph::vcount(graph)
    if ("name" %in% vertex_attr_names(graph))
        Z_names <- paste0(vertex_attr(graph, name="name"), "_d",
                          rep(1:d, each=nodes))
    else
        Z_names <- paste0("v", 1:nodes, "_d", rep(1:d, each=nodes))

    model <- structure(list(edges=edges, period=1, k=ncol(X), d=d,
                            ref=ref,
                            family=family, X=X,
                            beta_idx=1:ncol(X), beta_names=colnames(X),
                            Z_names=Z_names,
                            dropped=which(deg == 0), seed=seed,
                            verbose=verbose),
                       class="lsmfit")

    if (method == "MLE") {
        est <- lsm_MLE(theta, model, control=control)
    } else if (method == "MH") {
        theta <- lsm_MLE(theta, model, control=control)$par
        est   <- lsm_MH(theta, model, control=control)

        ## remove fixed values when ref units are used, apply names
        if (!is.null(ref)) {
            rm_idx <- which(apply(est$samples, 2, var) == 0)
            est$samples <- coda::mcmc(est$samples[,-rm_idx])
            varnames(est$samples) <- c(model$beta_names,
                                       model$Z_names[-rm_idx])
        } else {
            est$samples <- coda::mcmc(est$samples)
            varnames(est$samples) <- c(model$beta_names, model$Z_names)
        }
    } else {
        est <- NULL
    }

    model$method <- method
    model$graph <- graph
    model$estimate <- est

    model
}

build_model_data <- function(ft, graph)
{
    vars <- all.vars(ft)
    spec <- attr(ft, "specials")
    edges <- extract_edges(graph)

    if (!all(is.null(unlist(spec)))) {
        x <- lapply(names(spec), function(s) {
            if (!is.null(spec[[s]])) {
                fn <- get(s)
                fn(graph, vars[spec[[s]]])
            }
        })

        x <- do.call(cbind, x)
        x <- cbind(1, x)
    } else {
        x <- matrix(1, ncol=1, nrow=length(edges))
    }
    colnames(x) <- c("(Intercept)", attr(ft, "term.labels"))

    list(edges=edges, x=x)
}

## TODO: control should be used here
lsm_MLE <- function(theta, model, control)
{
    lwr <- c(rep(-Inf, ncol(model$X)), rep(-20, length(theta)-ncol(model$X)))
    upr <- c(rep(Inf, ncol(model$X)), rep(20, length(theta)-ncol(model$X)))
    est <- optim(theta, calc_likelihood, model=model, method="L-BFGS-B",
                 lower=lwr, upper=upr,
                 control=list(trace=model$verbose, fnscale=-1, maxit=500),
                 hessian=TRUE)

    est
}

lsm_MH <- function(theta, model, control=control.lsm(MCMC.burnin=2^10,
                                                     MCMC.samplesize=2^10,
                                                     MCMC.interval=10))
{
    Zstar <- matrix(theta[-model$beta_idx], ncol=model$d)
    if (!is.null(model$ref)) {
        Z     <- insert_ref(Zstar, model$ref, model$d)
        Z_idx <- (1:nrow(Z))[-model$ref$idx]
    } else {
        Z <- Zstar
        Z_idx <- 1:nrow(Zstar)
    }

    beta <- theta[model$beta_idx]

    .C_lsm_MH(model$edges, model$X, Z_idx, model$k, model$d,
              control$MCMC.burnin, control$MCMC.samplesize,
              control$MCMC.interval, beta, Z,
              model$family)
}

insert_ref <- function(pos, ref, d)
{
    n <- nrow(pos) + length(ref$idx)

    .C_insert_ref(ref$idx, ref$pos, (1:n)[-ref$idx], pos)
}

calc_likelihood <- function(theta, model=NULL)
{
    beta <- theta[model$beta_idx]
    X <- model$X

    Xbeta <- X %*% matrix(beta, ncol=1)

    if (!is.null(model$ref)) {
        pos <- insert_ref(matrix(theta[-model$beta_idx], ncol=model$d),
                          model$ref, model$d)
    } else {
        pos <- matrix(theta[-model$beta_idx], ncol=model$d)
    }

    D <- Xbeta - .C_dist_euclidean(pos)
    lik <- .C_llik_logit(model$edges, D)
    lik
}

start_random <- function(graph, ref, d)
{
    n <- igraph::vcount(graph)
    pos <- matrix(runif(n*d, min=-0.5, max=0.5), ncol=d, nrow=n)

    if (!is.null(ref))
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
           stop("unknown model family"))
}

procrustes <- function(Z, Zstar)
{
    d <- ncol(Z)
    Z <- scale(Z, center=TRUE, scale=FALSE)
    Ztran <- colMeans(Zstar)
    for (i in 1:d)
        Z[,i] <- Z[,i] + Ztran[i]

    A <- t(Z) %*% (Zstar %*% t(Zstar)) %*% Z
    E <- eigen(A, symmetric=TRUE)
    H <- E$vec %*% diag(sqrt(E$val)) %*% t(E$vec)

    t(t(Zstar) %*% Z %*% solve(H) %*% t(Z))
}
