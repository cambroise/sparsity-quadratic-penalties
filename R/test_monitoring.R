rm(list=ls())
library(quadrupen)
source("functions/simulations.R")

n <- 400
p <- 200
lambda2 <- 2

beta <- rep(rep(c(2,-2,0),c(15,15,70)),2)
data <- rlm(n, beta, corr=0.25, mu=3, r2=0.8)

subset <- sample(1:n, p/2)

fenchel <- elastic.net(data$x[subset, ], data$y[subset], lambda2=lambda2, nlambda1=6, control=list(monitor=2))

grandvt <- elastic.net(data$x[subset, ], data$y[subset], lambda2=lambda2, nlambda1=6, control=list(monitor=1))

delta_fenchel <- fenchel@monitoring$dist.to.opt
delta_fenchel [delta_fenchel < .Machine$double.eps] <- .Machine$double.eps

delta_grandvt <- grandvt@monitoring$dist.to.opt
delta_grandvt [delta_grandvt < .Machine$double.eps] <- .Machine$double.eps

delta_optimal <- grandvt@monitoring$dist.to.str
delta_optimal[delta_optimal < .Machine$double.eps] <- .Machine$double.eps

## plot(log10(delta_optimal), xlab="# iteration", type="l", lty=3, col="black", main="",
##      xlim=c(0,length(delta_optimal)), ylim=c(-8, max(log10(c(delta_optimal,delta_fenchel,delta_grandvt)))))
## lines(log10(delta_fenchel), xlab="# iteration", lty=2, col="red")
## lines(log10(delta_grandvt), xlab="# iteration", lty=1, col="blue")

pdf(file="bounds.pdf", width=12,height=7)
par(mar=c(2,2,0.1,0.1))
plot(log10(delta_optimal), ylab="", xlab="", type="l", lty=1, col="black", main="", lwd=2,
     xlim=c(0,length(delta_optimal)), ylim=c(-8, max(log10(c(delta_optimal,delta_fenchel,delta_grandvt)))))
lines(log10(delta_fenchel), xlab="# iteration", lty=3, col="red", lwd=2)
lines(log10(delta_grandvt), xlab="# iteration", lty=2, col="blue", lwd=2)


dev.off()
