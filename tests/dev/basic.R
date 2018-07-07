## =============================================================================
## Testing during development
## =============================================================================

options(width=120)

library(devtools)
library(lattice)
library(coda)
library(igraph)
setwd("~/Dropbox/lib/R/dynnet")
document("~/Dropbox/lib/R/dynnet")
load_all("~/Dropbox/lib/R/dynnet")

data(florentine, package="ergm")
(Flo  <- to_dynnet(flomarriage))

library(latentnet)

## -----------------------------------------------------------------------------
## Florentine without references
## -----------------------------------------------------------------------------

## tweak for better mixing
ctl <- control.lsm(MCMC.burnin=2^16, MCMC.interval=100)
model0 <- lsm(Flo ~ 1, ref=NULL, d=2, seed=1234)
system.time(model1 <- lsm(Flo ~ 1, ref=NULL, d=2, seed=1234, method="MH",
                          control=ctl))

ctl1 <- control.lsm(MCMC.burnin=2^16, MCMC.interval=100, dist_metric="euclidean")
model1 <- lsm(Flo ~ 1, ref=NULL, d=2, seed=1234, method="MH", control=ctl1)

ctl2 <- control.lsm(MCMC.burnin=2^16, MCMC.interval=100, dist_metric="euclidean2")
model2 <- lsm(Flo ~ 1, ref=NULL, d=2, seed=1234, method="MH", control=ctl2)

pdf("~/tmp/euclidean-test.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot_samples(model1, nsamp=1000, transformed=TRUE, cex=0.2, main="d")
plot_samples(model2, nsamp=1000, transformed=TRUE, cex=0.2, main="d2")
dev.off()


## -----------------------------------------------------------------------------
## Florentine With reference units
## -----------------------------------------------------------------------------

ref <- list(pos=matrix(c(-0.20, -0.20,
                          0.00,  0.20,
                          0.20, -0.20), ncol=2, byrow=TRUE),
            idx=c(7, 9, 14))

model2 <- lsm(Flo ~ 1, ref=ref, d=2, seed=1234)


document("~/lib/R/dynnet")
model3 <- lsm(Flo ~ 1, ref=ref, d=2, seed=1234, method="MH", control=ctl1)


ctl2 <- control.lsm(MCMC.burnin=2^20, MCMC.interval=200, MCMC.samplesize=2^12,
                    dist_metric="euclidean2")

mod <- lsm(Strike ~ 1, d=2, seed=1234, method="MH", control=ctl2)
mod.cov <- lsm(Strike ~ 1 + covmatch(group), d=2, seed=1234, method="MH", control=ctl2)

par(mfrow=c(1,2))
plot(mod)
plot(mod.cov)


plot_samples(mod.cov, nsamp=500, transformed=TRUE, cex=0.1)


par(mfrow=c(1,2))
plot(Strike)
plot(mod)


plot_mcmc(mod)
plot_mcmc(mod, transformed=TRUE)

par(mfrow=c(1,2))
plot(mod)
plot_samples(mod, nsamp=500, transformed=TRUE, cex=0.1)

X <- matrix(rnorm(200, sd=1), ncol=2)

## -----------------------------------------------------------------------------
## Strike data
## -----------------------------------------------------------------------------

data(Strike)

model0 <- lsm(Strike ~ 1, ref=NULL, d=2, seed=4321, verbose=TRUE)
summary(model0)

ctl <- control.lsm(MCMC.burnin=10^5, MCMC.interval=100, MCMC.samplesize=10^4,
                  dist_metric="euclidean2")
model1 <- lsm(Strike ~ 1, d=2, seed=4321, method="MH", control=ctl)
summary.lsmfit(model1)
plot.lsmfit(model1)
plot_mcmc(model1)
plot_mcmc(model1)

?predict.lsmfit(model1)
library(dynnet)
