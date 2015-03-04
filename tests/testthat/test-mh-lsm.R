context{"Metropolis-Hastings Latent Space Model"}

n <- 40
N <- n * (n-1) / 2
X <- matrix(1, nrow=N, ncol=1)
b <- c(0.25)
ref <- matrix(c(-1.0,0,0,1.0,1.0,0), ncol=2, nrow=3, byrow=TRUE)
S <- simulate_network(n, k=2, groups=3, periods=1, X=X, beta=b,
                      ref_pos=ref)
## mean_pos=ref,
## sigma_pos=rep(list(diag(2)*0.25), 3)

ref_idx <- S$DGP$ref_idx
z0 <- positions(S)
z0[-ref_idx,] <- 0.0
fam <- S$DGP$family
ITER <- 500000

## Rprof("~/tmp/mh.out", line.profiling=TRUE)
samp <- lsm(S$SIM, 0, z0, ref_idx, fam, iter=ITER)
## Rprof(NULL)

## summaryRprof("~/tmp/mh.out", lines="show")

par(mfrow=c(2,2))
plot(samp[[1]][,1], type="l")
plot(samp[[1]][,5], type="l")
plot(samp[[1]][,6], type="l")
plot(samp[[1]][,7], type="l")

for (i in 1:ncol(samp[[1]])) {
    plot(samp[[1]][,i], type="l", xlab=i, yla="value")
    Sys.sleep(3)
}


plot(samp[[1]][,ncol(samp[[1]])], type="l")


keep <- seq(ITER-2000, ITER, by=10)
hist(samp[[1]][keep, 1])
abline(v=mean(samp[[1]][keep, 1]), col="blue")
(est_b <- mean(samp[[1]][keep, 1]))

z_idx <- 2:(ncol(samp[[1]])-1)
est_Z <- matrix(colMeans(samp[[1]][keep, z_idx]), ncol=2)
sqrt(mean(((positions(S) - est_Z)[-ref_idx,])^2))


## pdf(file="~/tmp/truth-est.pdf", width=8, height=8)
truth <- positions(S)
plot(est_Z, xlim=c(-2.5, 2.5), ylim=c(-2.5, 2.5), pch=19, col="gray")
arrows(truth[-ref_idx,1], truth[-ref_idx,2],
       est_Z[-ref_idx,1], est_Z[-ref_idx,2], length=0.1)
points(truth, col="blue", pch=19)
points(ref, col="red", pch=1, cex=1.5)
## dev.off()


var(as.vector(truth[-ref_idx,]))
var(as.vector(est_Z[-ref_idx,]))

library(network)
plot(network(as.matrix(S$SIM$adj[[1]]), directed=FALSE))
