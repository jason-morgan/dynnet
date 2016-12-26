## ----------------------------------------------------------------------------
## This file is part of dynnet
##
## Copyright (C) 2016 Jason W. Morgan <jason.w.morgan@gmail.com>
##
## dynnet and is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program.  If not, see <http://www.gnu.org/licenses/>.
##
## ----------------------------------------------------------------------------


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
