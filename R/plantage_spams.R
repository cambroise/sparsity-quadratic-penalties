rm(list=ls())
library(spams)
library(mvtnorm)
library(lars)

## seed qui fait planter spam : 1539, 7719
## seed <- 1539
## seed ok pour SPAMs : 1819
## seed <- 1819
seed <- sample(1:10000,1)
set.seed(seed)
cat("\nseed #",seed)

## SIMULATION SETTINGS
n <- 100
p <- 400
rho <- 0.5
mu  <- 3
cat("\nCorrelation setup: rho =", rho)
ns <- floor(0.25 * min(n,p))
w <- rep(0,p)
w[sample(1:p,ns)] <- sample(c(-1,1),ns,replace=TRUE)*runif(ns,1,2)

Sigma  <- matrix(rho, p,p)
diag(Sigma) <- 1
x <- rmvnorm(n, sigma=Sigma, method="chol")
## choose sigma to control the r2
r2 <- 0.75
sigma  <- sqrt((1-r2)/r2 * t(w) %*% Sigma %*% w)
epsilon <- rnorm(n) * sigma
y <- mu + x %*% w + epsilon

## Reference path : lars package
cat(" lars package...\n")
out.lars  <- lars(x,y, max.steps=min(p,n))
coef.lars <- out.lars$beta[1:min(p,n),]

cat(" spams (Lars)...\n")
gc()
ybar <- mean(y)
xbar <- colMeans(x)
normx <- sqrt(drop(colSums(scale(x,xbar,FALSE)^2)))
mu <- ybar
coef.spm <- t(spams.lasso(scale(y,ybar,FALSE), scale(x,xbar,normx),
                          return_reg_path=TRUE, lambda1=0,
                          max_length_path= min(n,p),
                          numThreads = 1, cholesky=FALSE)[[2]])
coef.spm  <- scale(coef.spm,FALSE,normx)

cat(" spams (Fista)...\n")
gc()
coef.spm.f <- c()
beta0 <- matrix(rep(0,p),p,1)
for (l in out.lars$lambda) {
  beta0 <- spams.fistaFlat(scale(y,ybar,FALSE), scale(x,xbar,normx), beta0,
                             loss="square", regul="l1", lambda1=l, tol=1e-12)
  coef.spm.f <- cbind(coef.spm.f,beta0)
}
coef.spm.f  <- scale(t(coef.spm.f),FALSE,normx)

par(mfrow=c(1,3))
matplot(log10(out.lars$lambda),coef.spm , type="l", main="SPAMs (LARS)")
matplot(log10(out.lars$lambda),coef.spm.f, type="l", main="SPAMs (FISTA)")
matplot(log10(out.lars$lambda),coef.lars, type="l", main="LARS")
