rm(list=ls())
library(quadrupen)
library(glmnet)
library(lars)
library(spams)
source("functions/simulations.R")

seed <- sample(1:10000,1)
## seed <- 1819
## seed qui fait planter spam : 2473
## seed qui fonctionne pour spam: 7755, floor(0.10 * min(n,p))
## seed qui fonctionne pour spam: 1819, floor(0.25 * min(n,p))
set.seed(seed)
cat("\nseed #",seed)

nsim <- 30

scenario <- "high.dim.small"

## High dimensional settings for Lasso run
if (scenario == "high.dim.small") {
  n <- 40
  p <- 100
}

## High dimensional settings for Lasso run
if (scenario == "high.dim.med") {
  n <- 200
  p <- 2000
}
## High dimensional settings for Lasso run
if (scenario == "high.dim.hig") {
  n <- 1000
  p <- 20000
}

## SIMULATION SETTINGS
rho <- 0.4
cat("\nSmall correlation setup: rho =", rho)
ns <- floor(0.25 * min(n,p))
w <- rep(0,p)
w[sample(1:p,ns)] <- sample(c(-1,1),ns,replace=TRUE)*runif(ns,1,2)
data <- rlm(n, w, corr=rho, mu=3, r2=0.75)

## Reference path : lars package
time.lars <- system.time(fit.lars <- lar2fit(lars(data$x,data$y, max.steps=min(p,n)), colMeans(data$x)))[3]
objc.lars <- objective(fit.lars,data$x,data$y)
lambda    <- fit.lars@lambda1

cat(" quadrupen...")
gc()
fit.quad <- lasso(data$x,data$y, lambda, control=list(timer=TRUE))
time.quad <- fit.quad@monitoring$external.timer
objc.quad <- objective(fit.quad,data$x,data$y)
prec.quad <- sqrt(sum( (objc.quad - objc.lars)^2))/length(lambda)
                  
cat(" spams (LARS)...")
gc()
fit.spam <- lasso.spam(data$x,data$y,lambda)
time.spam <- fit.spam@monitoring$external.timer
objc.spam <- objective(fit.spam,data$x,data$y)
prec.spam <- sqrt(sum( (objc.spam - objc.lars)^2))/length(lambda)

cat(" spams (FISTA)...\n")
gc()
fit.spamf <- lasso.spam(data$x,data$y,lambda,type="fista",tol=1e-6)
time.spamf <- fit.spamf@monitoring$external.timer
objc.spamf <- objective(fit.spamf,data$x,data$y)
prec.spamf <- sqrt(sum( (objc.spamf - objc.lars)^2))/length(lambda)

res <- rbind(c(time.lars, time.quad, time.spam, time.spamf),
              c(0, prec.quad, prec.spam, prec.spamf))

colnames(res) <- c("lars","quadrupen","spams", "prox")
rownames(res) <- c("time", "precision")
print(res)

## effet du nombre d'itération d'optimisation sur le temps de calcul
## plot(diff(fit.quad@monitoring$pensteps.timer),
##     col=c(1,diff(cumsum(fit.quad@monitoring$it.optim)[cumsum(fit.quad@monitoring$it.active)[-1]])))
