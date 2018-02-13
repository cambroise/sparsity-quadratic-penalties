source("header.R")
header("~/svn/quadrupen")
library(lars)
library(spams)
source("functions/simulations.R")

seed <- sample(1:10000,1)
set.seed(seed)
nsim <- 10

##scenario <- "low.dim"
scenario <- "high.dim"

## Low dimensional settings for Lasso run
if (scenario == "low.dim") {
  n <- 400
  p <- 200
  file.out <- "output/glmnetVSquadra_p=200,n=400"
}
## High dimensional settings for Lasso run
if (scenario == "high.dim") {
  n <- 100
  p <- 400
  file.out <- "output/glmnetVSquadra_p=400,n=100"
}
ns <- floor(0.25 * min(n,p))

## SMALL CORRELATION SETTINGS
rho <- 0.5
cat("\nSmall correlation setup: rho =", rho)
data <- list()
fit.lars <- list()
time.lars.low <- rep(NA,nsim)
for (s in 1:nsim) {
  w <- rep(0,p)
  w[sample(1:p,ns)] <- sample(c(-1,1),ns,replace=TRUE)*runif(ns,1,2)
  data[[s]] <- rlm(n, w, corr=rho, mu=3, r2=0.75)
  time.lars.low[s] <- system.time(fit.lars[[s]] <- lar2fit(lars(data[[s]]$x,data[[s]]$y,
                                                                max.steps=min(p,n)),
                                                           colMeans(data[[s]]$x)))[3]
}

cat(" quadra...")
time.craf.low <- rep(NA,nsim)
prec.craf.low <- rep(NA,nsim)
for (s in 1:nsim) {
  out <- test.lasso(data[[s]], 1e-8, fit.lars[[s]], "quadra")
  time.craf.low[s] <- out[1]
  prec.craf.low[s] <- out[2]
}

cat(" spams...")
time.spam.low <- rep(NA,nsim)
prec.spam.low <- rep(NA,nsim)
for (s in 1:nsim) {
  out <- test.lasso(data[[s]], 1e-8, fit.lars[[s]], "spams")
  time.spam.low[s] <- out[1]
  prec.spam.low[s] <- out[2]
}

## MEDIUM CORRELATION SETTINGS
rho <- 0.4
cat("\nMedium correlation setup: rho =", rho)
data <- list()
time.lars.med <- rep(NA,nsim)
fit.lars <- list()
for (s in 1:nsim) {
  w <- rep(0,p)
  w[sample(1:p,ns)] <- sample(c(-1,1),ns,replace=TRUE)*runif(ns,1,2)
  data[[s]] <- rlm(n, w, corr=rho, mu=3, r2=0.75)
  time.lars.med[s] <- system.time(fit.lars[[s]] <- lar2fit(lars(data[[s]]$x,data[[s]]$y,
                                                                max.steps=min(p,n)),
                                                           colMeans(data[[s]]$x)))[3]
}

cat(" quadra...")
time.craf.med <- rep(NA,nsim)
prec.craf.med <- rep(NA,nsim)
for (s in 1:nsim) {
  out <- test.lasso(data[[s]], 1e-8, fit.lars[[s]], "quadra")
  time.craf.med[s] <- out[1]
  prec.craf.med[s] <- out[2]
}

cat(" spams...")
time.spam.med <- rep(NA,nsim)
prec.spam.med <- rep(NA,nsim)
for (s in 1:nsim) {
  out <- test.lasso(data[[s]], 1e-8, fit.lars[[s]], "spams")
  time.spam.med[s] <- out[1]
  prec.spam.med[s] <- out[2]
}

## HIGH CORRELATION SETTINGS
rho <- 0.8
cat("\nHigh correlation setup: rho =", rho)
data <- list()
fit.lars <- list()
time.lars.high <- rep(NA,nsim)
for (s in 1:nsim) {
  w <- rep(0,p)
  w[sample(1:p,ns)] <- sample(c(-1,1),ns,replace=TRUE)*runif(ns,1,2)
  data[[s]] <- rlm(n, w, corr=rho, mu=3, r2=0.75)
  time.lars.high[s] <- system.time(fit.lars[[s]] <- lar2fit(lars(data[[s]]$x,data[[s]]$y,
                                                                max.steps=min(p,n)),
                                                           colMeans(data[[s]]$x)))[3]
}

cat(" quadra...")
time.craf.high <- rep(NA,nsim)
prec.craf.high <- rep(NA,nsim)
for (s in 1:nsim) {
  out <- test.lasso(data[[s]], 1e-8, fit.lars[[s]], "quadra")
  time.craf.high[s] <- out[1]
  prec.craf.high[s] <- out[2]
}

cat(" spams...")
time.spam.high <- rep(NA,nsim)
prec.spam.high <- rep(NA,nsim)
for (s in 1:nsim) {
  out <- test.lasso(data[[s]], 1e-8, fit.lars[[s]], "spams")
  time.spam.high[s] <- out[1]
  prec.spam.high[s] <- out[2]
}

## OUTPUTS
rm(data)
save.image(file=paste(file.out,".RData",sep=""))

