loglik_logit <- function(y, Xb, d)
{
    sum(dbinom(y, 1, plogis(Xb - d), log=TRUE))
}

loglik_poisson <- function(y, Xb, d)
{
    sum(dpois(y, exp(Xb - d), log=TRUE))
}
