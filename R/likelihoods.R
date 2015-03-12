llik_fn <- function(family)
{
    switch(family,
           "logit"   = llik_logit,
           "poisson" = llik_poisson,
           stop("unknown family"))
}

llik_logit <- function(y, lp)
{
    sum(dbinom(y, 1, plogis(lp), log=TRUE))
}

llik_poisson <- function(y, lp)
{
    sum(dpois(y, 1, exp(lp), log=TRUE))
}
