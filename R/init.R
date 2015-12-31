## ##' Initialize model coeffients and latent positions for latent space and path
## ##' models.
## ##'
## ##' Initialize model coeffients and latent positions for latent space and path
## ##' models.
## ##' @title Initialize Model Coeffients and Latent Positions
## ##' @param Y igraph object.
## ##' @param X Model matrix for exogenous covariates.
## ##' @param ref Specified reference positions.
## ##' @return Vector of starting values.
## ##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
## lsm_init_coef <- function(Y, X, ref)
## {
##     Z <- lsm_init_Z(Y, ref$pos, ref$idx)
##     b <- lsm_init_b(Y, X, Z)
##     c(b, Z)
## }

## lsm_init_Z <- function(Y, ref_pos, ref_idx)
## {
##     k <- ncol(ref_pos)

##     start <- matrix(rnorm(k*nrow(Y[])), ncol=k, nrow=nrow(Y[]))
##     start[ref_idx,] <- ref_pos
##     as.vector(start)
## }

## ## lsm_init_Z <- function(Y, ref_pos, ref_idx)
## ## {
## ##     k <- ncol(ref_pos)

## ##     ## Geodesic distances
## ##     D <- igraph::shortest.paths(Y)
## ##     ## D <- scale(D, center=TRUE, scale=FALSE)
## ##     ## diag(D) <- 0

## ##     ## MDS to positions
## ##     Z0 <- cmdscale(D, k)

## ##     T <- MCMCpack::procrustes(Z0[ref_idx,], ref_pos,
## ##                               translation=TRUE, dilation=TRUE)
## ##     Z.star <- T$s * Z0 %*% T$R

## ##     for (i in 1:ncol(Z.star))
## ##         Z.star[,i] <- Z.star[,i] + T$t[i,]

## ##     Z.star
## ## }

## lsm_init_b <- function(Y, X, Z)
## {
##     rep(0, ncol(X))
## }
