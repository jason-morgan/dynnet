context{"Metropolis-Hastings Latent Space Model"}

ff <- read.table(text="
0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0
0 0 0 0 0 1 1 0 1 0 0 0 0 0 0 0
0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0
0 0 0 0 0 0 1 0 0 0 1 0 0 0 1 0
0 0 1 0 0 0 0 0 0 0 1 0 0 0 1 0
0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 1 0 1 0 0 0 1 0 0 0 0 0 0 0 1
0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
1 1 1 0 0 0 0 0 0 0 0 0 1 1 0 1
0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
0 0 0 1 1 0 0 0 0 0 0 0 0 0 1 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 1
0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
0 0 0 1 1 0 0 0 0 0 1 0 1 0 0 0
0 0 0 0 0 0 1 0 1 0 0 0 1 0 0 0")

ff <- as.matrix(ff)
nzero <- (colSums(ff) != 0)
ff <- ff[nzero, nzero]
FF <- dynnet_adjacency(ff)

## Florentine Marriage
ref_idx <- c(4, 5, 9)
z0 <- matrix(0, nrow=FF$nodes, ncol=2)
z0[4,] <- c( 1.0, 0.0) * 4
z0[5,] <- c( 0.0,-1.0) * 4
z0[9,] <- c(-1.0, 0.0) * 4

## z0[4,] <- c( 8.745830, 1.5456906)
## z0[7,] <- c( 2.632981, 4.6399155)
## z0[9,] <- c(-3.689019, 0.2952693)


## z0 <- structure(c(-4.64381650450682, -4.27194974719969, -1.82116438886952,
## 8.74582998301636, 4.32514514628101, -8.98363853830007, 2.63298064476044,
## 4.9347158462205, -3.68901935384851, -14.0517744666704, 9.07813996411139,
## 2.31318329501667, -10.4334578664838, 6.83475152224597, 1.74559268595525,
## 6.99952541745413, -0.233328088261106, -6.97839948993591, 1.54569064178947,
## -8.59395915988472, -4.75805347875178, 4.63991554489835, 9.98416601230354,
## 0.295269273148275, 3.82671232945801, -4.62945429066789, -1.75387469247412,
## 2.578135464681, -3.70093618176661, 3.01102333611935), .Dim = c(15L, 2L))

plot_ties <- function(adj, est)
{
    for (j in 1:(nrow(adj)-1)) {
        for (i in (j+1):nrow(adj)) {
            if (adj[i,j] == 1) {
                x <- est[c(i,j),1]
                y <- est[c(i,j), 2]
                lines(x, y, col="black")
            }
        }
    }
}


load_all("dynnet")
est <- lsm_MLE(FF, k=2, family="logit")

(z1 <- scale(matrix(est$par[2:length(est$par)], ncol=2), scale=FALSE))
plot(z1, col="white")
plot_ties(ff, z1)
text(z1[,1], z1[,2], labels=1:15, col="blue")

est$par[1]
sqrt(sum(dist(z1)^2))


FF$X <- list(matrix(1, ncol=1, nrow=15*14/2))
fam <- "logit"
BURN <- 2e5
SAMP <- 1e5
samp <- lsm_MH(FF, est$par[1], z1, ref_idx, fam, burnin=BURN, nsamples=SAMP)

keep <- seq(1, ITER, by=100)
hist(samp[[1]][keep, 1])
abline(v=mean(samp[[1]][keep, 1]), col="blue")
(est_b <- mean(samp[[1]][keep, 1]))

z_idx <- 2:(ncol(samp[[1]])-1)
all <- matrix(samp[[1]][keep, z_idx], ncol=2)
est_Z <- matrix(colMeans(samp[[1]][keep, z_idx]), ncol=2)
ref <- z0[ref_idx,]

plot(all, cex=0.2, col="gray")
points(est_Z, pch=19, xlim=c(min(all[,1]), max(all[,1])),
       ylim=c(min(all[,2]), max(all[,2])))
points(est_Z, cex=0.2, col="gray")
plot_ties(ff, est_Z)
text(est_Z[,1], est_Z[,2], labels=1:15, pos=1)
points(est_Z[c(4,5,9),], col="red", pch=1, cex=1.5)
## points(z0, col="blue", cex=1)

sqrt(sum(dist(est_Z)^2))

par(mfrow=c(1,2))
plot(samp[[1]][,1], type="l", xlab="b", ylab="value")
plot(samp[[1]][,ncol(samp[[1]])], type="l", xlab="lp", ylab="value")


## -----------------------------------------------------------------------------
## Simulated network.

set.seed(1104)
n <- 40
N <- n * (n-1) / 2
X <- matrix(1, nrow=N, ncol=1)
b <- c(0.50)
ref <- matrix(c(-1.0,0,0,1.0,1.0,0), ncol=2, nrow=3, byrow=TRUE)
S <- simulate_network(n, k=2, periods=1, X=X, beta=b,
                      ref_pos=ref)

## MLE
load_all("dynnet")
set.seed(1212)
nw    <- network(as.matrix(S$SIM$adj[[1]]), directed=FALSE)
start <- network.layout.fruchtermanreingold(nw, layout.par=NULL)
start <- c(0, as.vector(start))

est <- lsm_MLE(S$SIM, k=2, start=start, family="logit")
z1 <- matrix(est$theta[2:length(est$theta)], ncol=2)

hist(dist(z1) - dist(positions(S)))
?geodist




(z1 <- scale(matrix(est$par[2:length(est$par)], ncol=2), scale=FALSE))
plot(z1, col="white")
plot_ties(S$SIM$adj[[1]], z1)
text(z1[,1], z1[,2], labels=1:15, col="blue")

est$par[1]
sqrt(sum(dist(z1)^2))

M <- S$SIM$adj[[1]]
S$SIM$X <- list(matrix(1, ncol=1, nrow=nrow(M)*(nrow(M)-1)/2))
z0 <- matrix(0)
BURN <- 2e5
SAMP <- 1e5
samp <- lsm_MH(S$SIM, 0, z1, ref_idx, fam, burnin=BURN, nsamples=SAMP)

keep <- seq(1, ITER, by=100)
hist(samp[[1]][keep, 1])
abline(v=mean(samp[[1]][keep, 1]), col="blue")
(est_b <- mean(samp[[1]][keep, 1]))

z_idx <- 2:(ncol(samp[[1]])-1)
all <- matrix(samp[[1]][keep, z_idx], ncol=2)
est_Z <- matrix(colMeans(samp[[1]][keep, z_idx]), ncol=2)
ref <- z0[ref_idx,]

plot(all, cex=0.2, col="gray")
points(est_Z, pch=19, xlim=c(min(all[,1]), max(all[,1])),
       ylim=c(min(all[,2]), max(all[,2])))
points(est_Z, cex=0.2, col="gray")
plot_ties(ff, est_Z)
text(est_Z[,1], est_Z[,2], labels=1:15, pos=1)
points(est_Z[c(4,5,9),], col="red", pch=1, cex=1.5)
## points(z0, col="blue", cex=1)



## -----------------------------------------------------------------------------

library(network)
library(latentnet)
plot(network(as.matrix(S$SIM$adj[[1]]), directed=FALSE))

lnmod <- ergmm(network(ff[nzero,nzero], directed=FALSE) ~ euclidean(2),
               verbose=TRUE)

plot(lnmod1$mcmc.mle$Z, col="white")
text(lnmod1$mcmc.mle$Z, labels=1:15)
text(lnmod$mcmc.mle$Z, labels=1:15, col="gray")
## points(lnmod1$mcmc.pmode$Z, col="blue", pch=19)
text(lnmod1$mcmc.pmode$Z, labels=1:15, col="blue")
## points(est_Z, col="red", pch=19)
text(est_Z, labels=1:15, col="red")
text(lnmod1$mkl$Z, labels=1:15, col="green")
