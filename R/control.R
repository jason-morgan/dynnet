##' Control parameters for latent space models.
##'
##' Control parameters for latent space models.
##' @title Latent Space Model Control Parameters
##' @param MCMC.burnin Integer. Burn-in for MCMC sampler. Defaults to 1024.
##' @param MCMC.interval Integer. Number of steps to be taken between saving
##'     samples. Defaults to 100.
##' @param MCMC.samplesize Integer. Number of samples from the posterior to be
##'     saved. Defaults to 1024.
##' @return List of control parameters.
##' @author Jason W. Morgan \email{jason.w.morgan@@gmail.com}
##' @export
control.lsm <- function(MCMC.burnin=2^10,
                        MCMC.interval=100,
                        MCMC.samplesize=2^10)
{
    sapply(names(formals(sys.function())), get,
           envir=sys.frame(sys.parent(0)),
           simplify=FALSE)
}
