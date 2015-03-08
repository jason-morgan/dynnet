##' Initialize model coeffients and latent positions for latent space and path
##' models.
##'
##' Initialize model coeffients and latent positions for latent space and path
##' models.
##' @title Initialize Model Coeffients and Latent Positions
##' @param Y Adjacency matrix
##' @param k Positive integer. Number of dimensions for the latent social space.
##' @param X Model matrix for exogenous covariates.
##' @return Vector of starting values.
##' @author Jason W. Morgan \email{jason.w.morgan@gmail.com}
init_coef <- function(Y, k, X, ref_pos)
{
    Z <- init_Z(Y, k)
    b <- init_b(Y, X, Z)
    c(b, as.vector(Z))
}

init_Z <- function(Y, k, ref_pos)
{
    ## Geodesic distances
    D <- sna::geodist(Y, inf.replace=nrow(Y), ignore.eval=FALSE)
    ## MDS to positions
    ## Rotate to references
    D
}

init_b <- function(Y, X, Z)
{
    rep(0, ncol(X))
}
