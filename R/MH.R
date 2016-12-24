## Moved to lsm.R
## lsm_MH <- function(theta, model, control=control.lsm(MCMC.burnin=2^10,
##                                                      MCMC.samplesize=2^10,
##                                                      MCMC.interval=10))
## {
##     alpha <- theta[1]

##     if (!is.null(model$ref)) {
##         Z     <- insert_ref(matrix(theta[-1], ncol=model$k), model$ref, model$k)
##         Z_idx <- (1:nrow(Z))[-model$ref$idx]
##     } else {
##         Z <- matrix(theta[-1], ncol=model$k)
##         Z_idx <- 1:nrow(Z)
##     }

##     X <- matrix(0.0)
##     beta <- c(0)

##     .C_lsm_MH(model$edges, X, Z_idx, model$k,
##               control$MCMC.burnin, control$MCMC.samplesize,
##               control$MCMC.interval, alpha, beta, Z,
##               model$family)
## }
