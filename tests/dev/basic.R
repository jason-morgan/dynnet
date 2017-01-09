## =============================================================================
## Testing during development
## =============================================================================

options(width=120)

library(devtools)
library(lattice)
library(coda)
library(igraph)
setwd("~/lib/R/dynnet")
document("~/lib/R/dynnet")
load_all("~/lib/R/dynnet")

data(florentine, package="ergm")
(Flo  <- to_dynnet(flomarriage))

library(latentnet)

## -----------------------------------------------------------------------------
## Without reference units
## -----------------------------------------------------------------------------

## tweak for better mixing
ctl <- control.lsm(MCMC.burnin=2^18, MCMC.interval=100)
model0 <- lsm(Flo ~ 1, ref=NULL, d=2, seed=1234)
system.time(model1 <- lsm(Flo ~ 1, ref=NULL, d=2, seed=1234, method="MH",
                          control=ctl))

## -----------------------------------------------------------------------------
## With reference units
## -----------------------------------------------------------------------------

ref <- list(pos=matrix(c(-0.50, -0.50,
                          0.00,  0.50,
                          0.50, -0.50), ncol=2, byrow=TRUE),
            idx=c(7, 9, 14))

model2 <- lsm(Flo ~ 1, ref=ref, d=2, seed=1234)
model3 <- lsm(Flo ~ 1, ref=ref, d=2, seed=1234, method="MH", control=ctl)
