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
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
init_coef <- function(Y, X, ref_pos, ref_idx)
{
    Z <- init_Z(Y, ref_pos, ref_idx)
    b <- init_b(Y, X, Z)
    c(b, as.vector(Z))
}

init_Z <- function(Y, ref_pos, ref_idx)
{
    k <- ncol(ref_pos)

    ## Geodesic distances
    D <- igraph::shortest.paths(Y)

    ## MDS to positions
    Z0 <- cmdscale(D, k)

    T <- MCMCpack::procrustes(Z0[ref_idx,], ref_pos,
                              translation=TRUE, dilation=TRUE)
    Z.star <- T$s * Z0 %*% T$R

    for (i in 1:ncol(Z.star))
        Z.star[,i] <- Z.star[,i] + T$t[i,]

    Z.star
}

init_b <- function(Y, X, Z)
{
    rep(0, ncol(X))
}
