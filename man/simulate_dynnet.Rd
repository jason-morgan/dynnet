% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{simulate_dynnet}
\alias{simulate_dynnet}
\title{Simulate Dynamic and Static Networks}
\usage{
simulate_dynnet(n, k = 1, periods = 1, vattr = NULL,
  ref_fn = default_ref_pos(k), seed = NULL, family = "bernoulli", ...)
}
\arguments{
\item{n}{Positive integer. Number of nodes in the network.}

\item{k}{Positive integer. Number of dimensions in the latent space.}

\item{periods}{Positive integer. Number of periods in the dynnamic
network. If \code{k=1}, a static network is generated.}

\item{vattr}{List of vertex attributes. One data.frame per period.}

\item{ref_fn}{Function to generate the reference positions. Function must
return a list reference positions and indices.}

\item{seed}{Set the seed before generating the networks to assure
replicability. Default is \code{NULL}.}

\item{family}{Character string in \code{c("bernoulli", "poisson")}.}

\item{...}{Further parameters to be passed to subsequent functions.}
}
\value{
A \code{dynsim} object.
}
\description{
Simulate dynammic and static networks.
}
\details{
Simulate dynammic and static networks.
}
\author{
Jason W. Morgan \email{jason.w.morgan@gmail.com}
}

