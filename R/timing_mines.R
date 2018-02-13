library(lars)
library(quadrupen)
library(parallel)
library(glmnet)
library(spams)
source("functions/simulations.R")
source("functions/plots.R")

seed <- sample(1:10000,1)
set.seed(seed)
mc.cores <- 4
scenario <- "high.dim.small"
## scenario <- "high.dim.med"
## scenario <- "high.dim.hig"

## High dimensional settings for Lasso run
if (scenario == "high.dim.small") {
  nsim <- 100
  n <- 40
  p <- 100
  thres <- 10^c(0,-0.5,-1,-2,-3)
}

## High dimensional settings for Lasso run
if (scenario == "high.dim.med") {
  nsim <- 30
  n <- 200
  p <- 1000
  thres <- 10^c(0,-0.5,-1,-2)
}

## High dimensional settings for Lasso run
if (scenario == "high.dim.hig") {
  nsim <- 10
  n <- 400
  p <- 10000
  ##  thres <- 10^c(0,-0.5,-1,-2)
  thres <- 10^c(0,-0.5,-1)
}

file.out <- paste("output/timing_mines_p=",p,",n=",n,"_nsim=",nsim,".RData",sep="")
ns <- floor(0.25 * min(n,p))

## FUNCTION FOR PARALLELIZATION
draw.one.dataset <- function(s, rho) {
  cat(" ",s)
  w <- rep(0,p)
  w[sample(1:p,ns)] <- sample(c(-1,1),ns,replace=TRUE)*runif(ns,1,2)
  return(rlm(n, w, corr=rho, mu=3, r2=0.75))
}

call.lars <- function(data) {
    cat("+")
    time <- system.time(fit <- lar2fit(lars(data$x,data$y,
                                            max.steps=min(dim(data$x)),
                                            use.Gram=ifelse(ncol(data$x)<500,TRUE,FALSE)),
                                       colMeans(data$x)))[3]
    return(list(time=time,fit=fit))
}


## SMALL CORRELATION SETTINGS
rho <- 0.1
cat("\n=================================")
cat("\nSmall correlation setup: rho =", rho)
data <- list()
cat("\nData generation...")
data <- mclapply(1:nsim, draw.one.dataset, mc.cores=mc.cores, rho=rho)
cat("\n\n")

cat("\nLars... ")
out <- mclapply(data, call.lars, mc.cores=mc.cores)
time.lars.low <- sapply(out, function(l) l$time)
fit.lars <- sapply(out, function(l) l$fit)
cat("\n\n")

call.test.lasso <- function(s, method, eps) {
  cat("",s)
  return(test.lasso(data[[s]], eps, fit.lars[[s]], method))
}

cat("\nquadra (LARS)... ")
out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="quadrupen", eps=1e-8, mc.cores=mc.cores))
time.craf.low <- out[,1]
prec.craf.low <- out[,2]
cat("\n\n")

cat("\nquadra (FISTA)... ")
time.prox.low <- matrix(NA,nrow=length(thres),ncol=nsim)
prec.prox.low <- matrix(NA,nrow=length(thres),ncol=nsim)
for (t in seq_along(thres)) {
    cat("\n\t Current threshold =",thres[t])
    out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="quadrupen.prox", eps=thres[t], mc.cores=mc.cores ))
    time.prox.low[t,] <- out[, 1]
    prec.prox.low[t,] <- out[, 2]    
}
cat("\n\n")

cat("\nquadra (CD)... ")
time.path.low <- matrix(NA,nrow=length(thres),ncol=nsim)
prec.path.low <- matrix(NA,nrow=length(thres),ncol=nsim)
for (t in seq_along(thres)) {
    cat("\n\t Current threshold =",thres[t])
    out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="quadrupen.path", eps=thres[t], mc.cores=mc.cores ))
    time.path.low[t,] <- out[, 1]
    prec.path.low[t,] <- out[, 2]    
}
cat("\n\n")


## MEDIUM CORRELATION SETTINGS
rho <- 0.4
cat("\n=================================")
cat("\nMedium correlation setup: rho =", rho)
data <- list()
cat("\nData generation...")
data <- mclapply(1:nsim, draw.one.dataset, mc.cores=mc.cores, rho=rho)
cat("\n\n")

cat("\nLars... ")
out <- mclapply(data, call.lars, mc.cores=mc.cores)
time.lars.med <- sapply(out, function(l) l$time)
fit.lars <- sapply(out, function(l) l$fit)
cat("\n\n")

cat("\nquadra (LARS)... ")
out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="quadrupen", eps=eps, mc.cores=mc.cores))
time.craf.med <- out[,1]
prec.craf.med <- out[,2]
cat("\n\n")

cat("\nquadra (FISTA)... ")
time.prox.med <- matrix(NA,nrow=length(thres),ncol=nsim)
prec.prox.med <- matrix(NA,nrow=length(thres),ncol=nsim)
for (t in seq_along(thres)) {
    cat("\n\t Current threshold =",thres[t])
    out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="quadrupen.prox", eps=thres[t], mc.cores=mc.cores ))
    time.prox.med[t,] <- out[, 1]
    prec.prox.med[t,] <- out[, 2]    
}
cat("\n\n")

cat("\nquadra (CD)... ")
time.path.med <- matrix(NA,nrow=length(thres),ncol=nsim)
prec.path.med <- matrix(NA,nrow=length(thres),ncol=nsim)
for (t in seq_along(thres)) {
    cat("\n\t Current threshold =",thres[t])
    out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="quadrupen.path", eps=thres[t], mc.cores=mc.cores ))
    time.path.med[t,] <- out[, 1]
    prec.path.med[t,] <- out[, 2]    
}
cat("\n\n")


## HIGH CORRELATION SETTINGS
rho <- 0.8
cat("\n=================================")
cat("\nHigh correlation setup: rho =", rho)
data <- list()
cat("\nData generation...")
data <- mclapply(1:nsim, draw.one.dataset, mc.cores=mc.cores, rho=rho)
cat("\n\n")

cat("\nLars... ")
out <- mclapply(data, call.lars, mc.cores=mc.cores)
time.lars.high <- sapply(out, function(l) l$time)
fit.lars <- sapply(out, function(l) l$fit)
cat("\n\n")

cat("\nquadra (LARS)... ")
out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="quadrupen", eps=eps, mc.cores=mc.cores))
time.craf.high <- out[,1]
prec.craf.high <- out[,2]
cat("\n\n")

cat("\nquadra (FISTA)... ")
time.prox.high <- matrix(NA,nrow=length(thres),ncol=nsim)
prec.prox.high <- matrix(NA,nrow=length(thres),ncol=nsim)
for (t in seq_along(thres)) {
    cat("\n\t Current threshold =",thres[t])
    out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="quadrupen.prox", eps=thres[t], mc.cores=mc.cores ))
    time.prox.high[t,] <- out[, 1]
    prec.prox.high[t,] <- out[, 2]    
}
cat("\n\n")

cat("\nquadra (CD)... ")
time.path.high <- matrix(NA,nrow=length(thres),ncol=nsim)
prec.path.high <- matrix(NA,nrow=length(thres),ncol=nsim)
for (t in seq_along(thres)) {
    cat("\n\t Current threshold =",thres[t])
    out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="quadrupen.path", eps=thres[t], mc.cores=mc.cores ))
    time.path.high[t,] <- out[, 1]
    prec.path.high[t,] <- out[, 2]    
}
cat("\n\n")

## OUTPUTS
rm(data)
                                        #
save.image(file=file.out)

pdf(file=paste(file.out,".pdf",sep=""))
plot.mines(file.out)
dev.off()
