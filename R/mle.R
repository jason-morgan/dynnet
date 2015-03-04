## Z_step <- function(Z, b=NULL, net=NULL, llik_fn=logit_loglik)
## {
##     cols <- to_columns(net, c(b, Z))
##     d <- calc_distances(calc_positions(net$stan$T,
##                                        cols$pos, cols$traj))
##     d <- do.call(c, d)
##     llik_fn(net$stan$y, b, d)
## }

## b_step <- function(b, Z=NULL, net=NULL, llik_fn=logit_loglik)
## {
##     cols <- to_columns(net, c(b, Z))    
##     d <- calc_distances(calc_positions(net$stan$T,
##                                        cols$pos, cols$traj))
##     d <- do.call(c, d)
##     llik_fn(net$stan$y, b, d)
## }

## data_lsm <- function(net)
## {
##     adj <- get_adjacency(net, period=1)    
## }

## est_mle_lsm <- function(net) {

##     data_lsm(net)

##     ## tol <- sqrt(.Machine$double.eps)
##     ## i   <- 1

##     ## current_b <- mean(net$stan$y)
##     ## current_Z <- rep(0, (net$stan$n-net$stan$K-1)*Net$stan$K * 2)

##     ## if (verbose > 0) {
##     ##     cat("## ", rep("-", 77), "\n", sep="")        
##     ##     cat("## Estimating starting positions and trajectories...\n")
##     ## }
            
##     ## est_start <- optim(current_Z, Z_step, b=current_b, net=Net, llik_fn=llik_fn,
##     ##                    method="SANN",
##     ##                    control=list(maxit=10000, fnscale=-1, trace=verbose-1))

##     ## current_Z <- est_start$par
##     ## current_llik <- est_start$value
##     ## old_lik <- current_llik + 10        # algorithm will run at least once
##     ## if (verbose > 0) cat("##   Initial log-likelihood:", current_llik, "\n")    
    
##     ## while (abs(current_llik - old_lik) > tol && i <= iter) {
##     ##     ## update Z, b fixed
##     ##     if (verbose > 0) {
##     ##         cat("## ", rep("-", 77), "\n", sep="")
##     ##         cat("## Iteration:", i, "\n##   Z step...\n")
##     ##     }
##     ##     est_Z <- optim(current_Z, Z_step, b=current_b, net=Net, llik_fn=llik_fn,
##     ##                    method=method,
##     ##                    control=list(maxit=10000, fnscale=-1, trace=verbose-1))

##     ##     current_Z <- est_Z$par

##     ##     ## update b, Z fixed
##     ##     if (verbose > 0) cat("##   b step...\n")
##     ##     est_b <- optim(current_b, b_step, Z=current_Z, net=Net, llik_fn=llik_fn,
##     ##                    method="BFGS",
##     ##                    control=list(maxit=10000, fnscale=-1, trace=verbose-1))

##     ##     current_b <- est_b$par

##     ##     old_llik     <- current_llik
##     ##     current_llik <- est_b$value
##     ##     i <- i+1
##     ##     if (verbose > 0) cat("##   New log-likelihood:", current_llik, "\n")
##     ## }

##     ## if (verbose > 0) {
##     ##     cat("## ", rep("-", 77), "\n", sep="")
##     ##     cat("##   Final log-likelihood:", current_llik, "\n")
##     ##     cat("##", rep("-", 77), "\n")        
##     ## }
    
##     ## list(b=current_b, Z=current_Z)
## }
