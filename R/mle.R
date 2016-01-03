
## lsm_MLE <- function(model, verbose=TRUE)
## {
##     model$llik <- mk_log_likelihood(model$family, k=model$k)

##     ## Hacky -- get the response
##     y <- get_adjacency(model$network, period=model$period)
##     y <- as.matrix(Matrix::tril(y))
##     model$y <- y[lower.tri(y)]

##     ## Set up indices for the variables to be estimated.
##     model$b_idx    <- 1:ncol(model$X)
##     model$Z_idx    <- (length(model$b_idx)+1):(length(model$start))
##     model$template <- matrix(model$start[model$Z_idx], ncol=model$k)

##     theta <- c(model$start[model$b_idx],
##                as.vector(matrix(model$start[model$Z_idx],
##                                 ncol=model$k)[-model$ref$idx,]))

##     print(theta)

##     model$Z_idx <- (length(model$b_idx)+1):(length(theta))

##     ## ## Initial estimates of Z
##     if (isTRUE(verbose)) cat("Calculating initial locations...\n")
##     cur <- optim(theta[model$Z_idx], MLE_Z_update, theta=theta, model=model,
##                  ## method="SANN",
##                  method="BFGS",
##                  control=list(trace=0, fnscale=-1, maxit=100000))
##     theta[model$Z_idx] <- cur$par

##     old_lik <- 0
##     new_lik <- 1
##     iter    <- 0
##     MAXITER <- 10
##     tol     <- sqrt(.Machine$double.eps)

##     while (MAXITER > iter && (abs(new_lik - old_lik) > tol)) {
##         old_lik <- new_lik

##         cur   <- MLE_est_1step(theta, model)
##         ## cur   <- MLE_est_2step(theta, model)
##         theta <- cur$theta

##         if (isTRUE(verbose) && iter %% 5 == 0) {
##             tmpd <- mean(dist(matrix(theta[model$Z_idx], ncol=model$k)))
##             cat("Iteration:", iter, "\n")
##             cat("  beta:", theta[model$b_idx], "\n")
##             cat("  mean distance:", tmpd, "\n")
##             cat("  log-likelihood:", cur$value, "\n")
##         }

##         iter    <- iter + 1
##         new_lik <- cur$est$value
##     }

##     list(theta=center_Z(theta, model$b_idx, model$Z_idx, model$k),
##          last=cur, model=model)

##     list(theta=theta, last=cur, model=model)
## }

## MLE_est_1step <- function(theta, model)
## {
##     cur <- optim(theta, model$llik, model=model,
##                  ## method="BFGS",
##                  control=list(trace=0, fnscale=-1))
##     theta <- cur$par
##     list(theta=theta, est=cur)
## }

## MLE_est_2step <- function(theta, model)
## {
##     cur <- optim(theta[model$b_idx], MLE_b_update, theta=theta, model=model,
##                  method="Brent", lower=-100, upper=100,
##                  control=list(trace=0, fnscale=-1))
##     theta[model$b_idx] <- cur$par

##     cur <- optim(theta[model$Z_idx], MLE_Z_update, theta=theta, model=model,
##                  ## method="BFGS",
##                  control=list(trace=0, fnscale=-1))
##     theta[model$Z_idx] <- cur$par

##     list(theta=theta, est=cur)
## }

## MLE_Z_update <- function(Z, theta=NULL, model=NULL)
## {
##     theta[model$Z_idx] <- Z
##     model$llik(theta, model)
## }

## MLE_b_update <- function(b, theta=NULL, model=NULL)
## {
##     theta[model$b_idx] <- b
##     model$llik(theta, model)
## }

## center_Z <- function(theta, b_idx, Z_idx, k)
## {
##     Z_est <- scale(matrix(theta[Z_idx], ncol=k), scale=FALSE)
##     c(theta[b_idx], as.vector(Z_est))
## }

## make_Z <- function(theta, Z_idx, template, ref_idx)
## {
##     template[-ref_idx,] <- matrix(theta[Z_idx], ncol=ncol(template))
##     template
## }

## mk_log_likelihood <- function(family, k)
## {
##     llik <- llik_fn(family)

##     function(theta, model, lambda=0.25)
##     {
##         b  <- theta[model$b_idx]
##         Z  <- make_Z(theta, model$Z_idx, model$template, model$ref$idx)
##         d  <- as_distance_vector(Z)
##         lp <- (model$X %*% b) - d

##         ## norm <- apply(Z, 1, function(r) sqrt(sum(r^2)))
##         ## norm <- norm_euclidean(Z)
##         ## penalty <- distance_penalty(Z)

##         l <- llik(model$y, lp)
##         p <- lambda * abs(b)

##         ## llik(model$y, lp) - lambda * sum(log(norm))
##         ## llik(model$y, lp) - lambda * penalty
##         l - p
##     }
## }
