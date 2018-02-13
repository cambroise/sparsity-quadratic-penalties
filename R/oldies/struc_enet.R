rm(list=ls())
gc()
library(elasticnet)
library(R.utils)
sourceDirectory("functions")

seed <- sample(1:10000,1)
set.seed(seed)

p <- 100
n <- 200
beta <- rep(c(0,2,0,-2),each=5)
beta <- rep(beta,p%/%20)

data <- example_c_enet(n, beta, cor=0.25)
x <- data$x
y <- data$y

D <- matrix(0,p,p)
for (k in 1:20) {
  D[((k-1)*5+1):(k*5),((k-1)*5+1):(k*5)] <- 0.9
}
diag(D) <- 1

out.lass <- crafter(x, y, gamma=0)
out.enet <- crafter(x, y, gamma=1)
out.stru <- crafter(x, y, gamma=1, D)

par(mfrow=c(1,3))
plot(out.stru$fit, main="Structured Enet, gamma=10",
     xvar="fraction", col=rep(1:20,each=5), lty=rep(1:20,each=5))
plot(out.enet$fit, main="Elastic-net, gamma=10",
     xvar="fraction", col=rep(1:20,each=5), lty=rep(1:20,each=5))
plot(out.lass$fit, main="Lasso",
     xvar="fraction", col=rep(1:20,each=5), lty=rep(1:20,each=5))
