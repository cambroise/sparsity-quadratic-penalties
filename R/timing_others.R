library(quadrupen)
library(lars)
library(glmnet)
library(spams)
source("functions/simulations.R")
source("functions/plots.R")

seed <- sample(1:10000,1)
set.seed(seed)

##scenario <- "high.dim.med"
##scenario <- "high.dim.hig"
scenario <- "high.dim.small"

## High dimensional settings for Lasso run
if (scenario == "high.dim.small") {
  nsim <- 10
  n <- 40
  p <- 100
}

## High dimensional settings for Lasso run
if (scenario == "high.dim.med") { 
  nsim <- 30
  n <- 200
  p <- 1000
}
## High dimensional settings for Lasso run
if (scenario == "high.dim.hig") {
  nsim <- 10
  n <- 400
  p <- 10000
}

file.out <- paste("output/timing_others_p=",p,",n=",n,"_nsim=",nsim,".RData",sep="")

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

cat("\nSPAMs (LARS)... ")
out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="spams.lars", eps=1e-8, mc.cores=mc.cores))
time.spam.low <- out[,1]
prec.spam.low <- out[,2]
cat("\n\n")

cat("\nglmnet... ")
thres.glmn <- 10^c(-5, -8, -10, -12, -14, -16, -18, -20)
time.glmn.low <- matrix(NA,nrow=length(thres.glmn),ncol=nsim)
prec.glmn.low <- matrix(NA,nrow=length(thres.glmn),ncol=nsim)
for (t in seq_along(thres.glmn)) {
    cat("\n\t Current threshold =",thres[t])
    out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="glmnet", eps=thres.glmn[t], mc.cores=mc.cores ))
    time.glmn.low[t,] <- out[, 1]
    prec.glmn.low[t,] <- out[, 2]    
}
cat("\n\n")

cat("\nSPAMs (fista)... ")
thres.spams <- 10^c(-1, -4, -5)
time.spamf.low <- matrix(NA,nrow=length(thres.spams),ncol=nsim)
prec.spamf.low <- matrix(NA,nrow=length(thres.spams),ncol=nsim)
for (t in seq_along(thres.spams)) {
    cat("\n\t Current threshold =",thres[t])
    out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="spams.prox", eps=thres.spams[t], mc.cores=mc.cores ))
    time.spamf.low[t,] <- out[, 1]
    prec.spamf.low[t,] <- out[, 2]    
}
cat("\n\n")

## MEDIUM CORRELATION SETTINGS
rho <- 0.4
cat("\n=================================")
cat("\nMedium correlation setup: rho =", rho)
cat("\nData generation...")
data <- mclapply(1:nsim, draw.one.dataset, mc.cores=mc.cores, rho=rho)
cat("\n\n")

cat("\nLars... ")
out <- mclapply(data, call.lars, mc.cores=mc.cores)
time.lars.med <- sapply(out, function(l) l$time)
fit.lars <- sapply(out, function(l) l$fit)
cat("\n\n")

call.test.lasso <- function(s, method, eps) {
  cat("",s)
  return(test.lasso(data[[s]], eps, fit.lars[[s]], method))
}

cat("\nquadra (LARS)... ")
out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="quadrupen", eps=1e-8, mc.cores=mc.cores))
time.craf.med <- out[,1]
prec.craf.med <- out[,2]
cat("\n\n")

cat("\nSPAMs (LARS)... ")
out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="spams.lars", eps=1e-8, mc.cores=mc.cores))
time.spam.med <- out[,1]
prec.spam.med <- out[,2]
cat("\n\n")

cat("\nglmnet... ")
thres.glmn <- 10^c(-5, -7, -9, -11, -13) 
time.glmn.med <- matrix(NA,nrow=length(thres.glmn),ncol=nsim)
prec.glmn.med <- matrix(NA,nrow=length(thres.glmn),ncol=nsim)
for (t in seq_along(thres.glmn)) {
    cat("\n\t Current threshold =",thres[t])
    out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="glmnet", eps=thres.glmn[t], mc.cores=mc.cores ))
    time.glmn.med[t,] <- out[, 1]
    prec.glmn.med[t,] <- out[, 2]    
}
cat("\n\n")

cat("\nSPAMs (fista)... ")
thres.spams <- 10^c(-1, -4, -5)
time.spamf.med <- matrix(NA,nrow=length(thres.spams),ncol=nsim)
prec.spamf.med <- matrix(NA,nrow=length(thres.spams),ncol=nsim)
for (t in seq_along(thres.spams)) {
    cat("\n\t Current threshold =",thres[t])
    out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="spams.prox", eps=thres.spams[t], mc.cores=mc.cores ))
    time.spamf.med[t,] <- out[, 1]
    prec.spamf.med[t,] <- out[, 2]    
}
cat("\n\n")

## HIGH CORRELATION SETTINGS
rho <- 0.8
cat("\n=================================")
cat("\nHigh correlation setup: rho =", rho)
cat("\nData generation...")
data <- mclapply(1:nsim, draw.one.dataset, mc.cores=mc.cores, rho=rho)
cat("\n\n")

cat("\nLars... ")
out <- mclapply(data, call.lars, mc.cores=mc.cores)
time.lars.high <- sapply(out, function(l) l$time)
fit.lars <- sapply(out, function(l) l$fit)
cat("\n\n")

call.test.lasso <- function(s, method, eps) {
  cat("",s)
  return(test.lasso(data[[s]], eps, fit.lars[[s]], method))
}

cat("\nquadra (LARS)... ")
out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="quadrupen", eps=1e-8, mc.cores=mc.cores))
time.craf.high <- out[,1]
prec.craf.high <- out[,2]
cat("\n\n")

cat("\nSPAMs (LARS)... ")
out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="spams.lars", eps=1e-8, mc.cores=mc.cores))
time.spam.high <- out[,1]
prec.spam.high <- out[,2]
cat("\n\n")

cat("\nglmnet... ")
thres.glmn <- 10^c(-5, -7, -8, -9, -10)
time.glmn.high <- matrix(NA,nrow=length(thres.glmn),ncol=nsim)
prec.glmn.high <- matrix(NA,nrow=length(thres.glmn),ncol=nsim)
for (t in seq_along(thres.glmn)) {
    cat("\n\t Current threshold =",thres[t])
    out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="glmnet", eps=thres.glmn[t], mc.cores=mc.cores ))
    time.glmn.high[t,] <- out[, 1]
    prec.glmn.high[t,] <- out[, 2]    
}
cat("\n\n")

cat("\nSPAMs (fista)... ")
thres.spams <- 10^c(-1, -2, -4)
time.spamf.high <- matrix(NA,nrow=length(thres.spams),ncol=nsim)
prec.spamf.high <- matrix(NA,nrow=length(thres.spams),ncol=nsim)
for (t in seq_along(thres.spams)) {
    cat("\n\t Current threshold =",thres[t])
    out <- do.call(rbind, mclapply(1:nsim, call.test.lasso, method="spams.prox", eps=thres.spams[t], mc.cores=mc.cores ))
    time.spamf.high[t,] <- out[, 1]
    prec.spamf.high[t,] <- out[, 2]    
}
cat("\n\n")

## OUTPUTS
rm(data)
save.image(file=file.out)

pdf(file=paste(file.out,".pdf",sep=""))
plot.others(file.out)
dev.off()
