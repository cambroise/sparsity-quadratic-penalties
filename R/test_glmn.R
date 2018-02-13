library(quadrupen)
library(lars)
library(glmnet)
source("functions/simulations.R")

seed <- sample(1:10000,1)
set.seed(seed)
nsim <- 5

thres1 <- 1e-3
thres2 <- 1e-5
thres3 <- 1e-7
thres4 <- 1e-9
thres5 <- 1e-12

scenario <- "low.dim"
##scenario <- "high.dim"

## Low dimensional settings for Lasso run
if (scenario == "low.dim") {
  n <- 800
  p <- 200
  file.out <- "output/glmnetVScrafter_p=200,n=400"
}
## High dimensional settings for Lasso run
if (scenario == "high.dim") {
  n <- 100
  p <- 400
  file.out <- "output/glmnetVScrafter_p=400,n=100"
}

rho <- 0.1
cat("\nCorrelation setup: rho =", rho)

data <- rlm(n, runif(p), corr=rho, r2=0.8, sparse=FALSE, prob=0.01, size=5)
y <- data$y
x <- data$x
xbar  <- colMeans(x)
normx <- sqrt(drop(colSums(scale(x,xbar,FALSE)^2)))
fit.lars <- lar2fit(lars(x,y), xbar)

lambda <- fit.lars@lambda1
fit.glmn1 <- glm2fit(glmnet(x,y,lambda=lambda/sqrt(n),thres=thres1),n)
fit.glmn2 <- glm2fit(glmnet(x,y,lambda=lambda/sqrt(n),thres=thres2),n)
fit.glmn3 <- glm2fit(glmnet(x,y,lambda=lambda/sqrt(n),thres=thres3),n)
fit.glmn4 <- glm2fit(glmnet(x,y,lambda=lambda/sqrt(n),thres=thres4),n)
fit.glmn5 <- glm2fit(glmnet(x,y,lambda=lambda/sqrt(n),thres=thres5),n)
fit.quad <- lasso(x,y,lambda1=lambda)

## remove last lambda from all ft (instable)
nlbd <- length(lambda)
fit.glmn1@coefficients <- fit.glmn1@coefficients[-nlbd, ]
fit.glmn2@coefficients <- fit.glmn2@coefficients[-nlbd, ]
fit.glmn3@coefficients <- fit.glmn3@coefficients[-nlbd, ]
fit.glmn4@coefficients <- fit.glmn4@coefficients[-nlbd, ]
fit.glmn5@coefficients <- fit.glmn5@coefficients[-nlbd, ]

fit.quad@coefficients <- fit.quad@coefficients[-nlbd, ]
fit.lars@coefficients <- fit.lars@coefficients[-nlbd, ]
fit.glmn1@lambda1 <- fit.glmn1@lambda1[-nlbd]
fit.glmn2@lambda1 <- fit.glmn2@lambda1[-nlbd]
fit.glmn3@lambda1 <- fit.glmn3@lambda1[-nlbd]
fit.glmn4@lambda1 <- fit.glmn4@lambda1[-nlbd]
fit.glmn5@lambda1 <- fit.glmn5@lambda1[-nlbd]
fit.quad@lambda1 <- fit.quad@lambda1[-nlbd]
fit.lars@lambda1 <- fit.lars@lambda1[-nlbd]
fit.glmn1@mu <- fit.glmn1@mu[-nlbd]
fit.glmn2@mu <- fit.glmn2@mu[-nlbd]
fit.glmn3@mu <- fit.glmn3@mu[-nlbd]
fit.glmn4@mu <- fit.glmn4@mu[-nlbd]
fit.glmn5@mu <- fit.glmn5@mu[-nlbd]
fit.quad@mu <- fit.quad@mu[-nlbd]
fit.lars@mu <- fit.lars@mu[-nlbd]

obj.lars <- objective(fit.lars, x,y, mu=data$mu)
obj.glmn1 <- objective(fit.glmn1, x,y, mu=data$mu)
obj.glmn2 <- objective(fit.glmn2, x,y, mu=data$mu)
obj.glmn3 <- objective(fit.glmn3, x,y, mu=data$mu)
obj.glmn4 <- objective(fit.glmn4, x,y, mu=data$mu)
obj.glmn5 <- objective(fit.glmn5, x,y, mu=data$mu)
obj.quad <- objective(fit.quad, x,y, mu=data$mu)

dif.glmn1 <- fit.glmn1 
dif.glmn1@coefficients <- abs(dif.glmn1@coefficients - fit.lars@coefficients)
dif.glmn2 <- fit.glmn2 
dif.glmn2@coefficients <- abs(dif.glmn2@coefficients - fit.lars@coefficients)
dif.glmn3 <- fit.glmn3 
dif.glmn3@coefficients <- abs(dif.glmn3@coefficients - fit.lars@coefficients)
dif.glmn4 <- fit.glmn4 
dif.glmn4@coefficients <- abs(dif.glmn4@coefficients - fit.lars@coefficients)
dif.glmn5 <- fit.glmn5 
dif.glmn5@coefficients <- abs(dif.glmn5@coefficients - fit.lars@coefficients)

dif.quad <- fit.quad
dif.quad@coefficients <- abs(dif.quad@coefficients - fit.lars@coefficients)

xvar <- log10(lambda)[-nlbd]

par(mfrow=c(2,2))
plot(xvar, obj.lars,type="l",col="black", xlab="log10(lambda)", ylab="fonction objectif", )
lines(xvar, obj.glmn1,col="blue",lty=1)
lines(xvar, obj.glmn2,col="blue",lty=2)
lines(xvar, obj.glmn3,col="blue",lty=3)
lines(xvar, obj.glmn4,col="blue",lty=4)
lines(xvar, obj.glmn5,col="blue",lty=5)
lines(xvar, obj.quad,col="red")
legend("topleft", legend=c("glmnet 1e-3","glmnet 1e-5","glmnet 1e-7","glmnet 1e-9","glmnet 1e-12",
                    "crafter"), lty=c(1:5,1), col=c(rep("blue",5), "red"),bty="n")

plot(xvar, log10(cumsum(abs(obj.glmn1-obj.lars))),
     xlab="log10(lambda)", ylab="Écart cumulé objectif (log10)", type="l",col="blue", lty=1)
lines(xvar, log10(cumsum(abs(obj.glmn2-obj.lars))),col="blue", lty=2)
lines(xvar, log10(cumsum(abs(obj.glmn3-obj.lars))),col="blue", lty=3)
lines(xvar, log10(cumsum(abs(obj.glmn4-obj.lars))),col="blue", lty=4)
lines(xvar, log10(cumsum(abs(obj.glmn5-obj.lars))),col="blue", lty=5)
lines(xvar, log10(cumsum(abs(obj.quad-obj.lars))),col="red")

plot(xvar,log10(cumsum(rowSums(abs(fit.glmn1@intercept -fit.lars@intercept) + (fit.glmn1@coefficients-fit.lars@coefficients)^2))),
     xlab="log10(lambda)", ylab="Écart cumulé coefficients (log10)", type="l",col="blue", lty=1)
lines(xvar, log10(cumsum(abs(fit.glmn2@intercept -fit.lars@intercept) + rowSums((fit.glmn2@coefficients-fit.lars@coefficients)^2))),col="blue", lty=2)
lines(xvar, log10(cumsum(abs(fit.glmn3@intercept -fit.lars@intercept) + rowSums((fit.glmn3@coefficients-fit.lars@coefficients)^2))),col="blue", lty=3)
lines(xvar, log10(cumsum(abs(fit.glmn4@intercept -fit.lars@intercept) + rowSums((fit.glmn4@coefficients-fit.lars@coefficients)^2))),col="blue", lty=4)
lines(xvar, log10(cumsum(abs(fit.glmn5@intercept -fit.lars@intercept) + rowSums((fit.glmn5@coefficients-fit.lars@coefficients)^2))),col="blue", lty=5)

lines(xvar, log10(cumsum(abs(fit.quad@intercept -fit.lars@intercept) + rowSums((fit.quad@coefficients-fit.lars@coefficients)^2))),col="red")

plot(xvar, log10(cumsum(colSums((predict(fit.glmn1,x) -predict(fit.lars,x))^2))),
     type="l",col="blue", lty=1, xlab="log10(lambda)", ylab="Écart cumulé prédiction log10")
lines(xvar, log10(cumsum(colSums((predict(fit.glmn2,x) -predict(fit.lars,x))^2))),col="blue", lty=2)
lines(xvar, log10(cumsum(colSums((predict(fit.glmn3,x) -predict(fit.lars,x))^2))),col="blue", lty=3)
lines(xvar, log10(cumsum(colSums((predict(fit.glmn4,x) -predict(fit.lars,x))^2))),col="blue", lty=4)
lines(xvar, log10(cumsum(colSums((predict(fit.glmn5,x) -predict(fit.lars,x))^2))),col="blue", lty=5)
lines(xvar, log10(cumsum(colSums((predict(fit.quad,x) -predict(fit.lars,x))^2))),col="red")
title(main="\n\nÉcart à la solution optimale (package lars)", outer=TRUE)
