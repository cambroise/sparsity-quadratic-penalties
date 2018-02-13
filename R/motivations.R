rm(list=ls())
gc()
source("functions/simulations.R")
source("functions/plots.R")

## COLOR BLIND PALETTE
library(ggplot2)
library(scales)
library(Hmisc)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scale_colour_manual(values=cbbPalette)

seed <- sample(1:10000,1)
set.seed(seed) ## 4157
## Sparsity settings
sp.hig <- sample(rep(c(1,-1,0),c(5,5,90)))
sp.med <- sample(rep(c(1,-1,0),c(15,15,70)))
sp.low <- sample(rep(c(1,-1,0),c(20,20,60)))
beta <- sp.med

p <- length(beta)
prop.test <- 9/10
r2  <- 0.8
thres.glmn <- 10^c(-9, -4, -1)
K <- 100
cor <- 0.8

test <- list()
supp <- list()
prec <- list()
time <- list()
out.all <- c()
l <- 0
lambda <- 10^seq(1, -.8, len=100)

for (n0 in c(p/2, p, 2*p)) {
  cat("\n\n========================================================================")
  cat("\nN = ",n0)
  cat("\n========================================================================")
  l <- l+1
  n <- n0*10
  
  ## METHODS COMPARISON
  out <- outsample.test.error(K, n, beta, cor, r2, lambda, thres.glmn, prop.test=prop.test)
  
  ## Graphical outputs
  method <- out$method
  method <- rep(NA, length(out$method))
  method[out$method == "quad"] <- "quadrupen"
  method[out$method == "glmn 1e-09"] <- "glmnet high"
  method[out$method == "glmn 1e-04"] <- "glmnet med"
  method[out$method == "glmn 0.1"] <- "glmnet low"
  method <- factor(method, level=c("quadrupen", "glmnet low", "glmnet med", "glmnet high"), ordered=TRUE)
  out$method <- method
  out$np.ratio <- rep(paste("n/p =",n0/p), length(nrow(out)))

  out.all <- rbind(out.all, out)

  ## test[[l]] <- 
  ##   ggplot(out, aes(x=lambda,y=sqrt(error.mse),colour=method,group=method,lty=method)) +
  ##     stat_summary(fun.data="mean_cl_normal", geom="smooth", width=3, alpha=0.1) +
  ##       scale_x_continuous(trans=log2_trans())
  
  ## supp[[l]] <- 
  ##   ggplot(out, aes(x=lambda, y=error.sup, colour=method, group=method,lty=method)) +
  ##     stat_summary(fun.data="mean_cl_normal", geom="smooth", width=3, alpha=0.1)+
  ##       scale_x_continuous(trans=log2_trans())
}

test.all <- 
  ggplot(out.all, aes(x=lambda,y=sqrt(error.mse),colour=method,group=method,lty=method)) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth", width=3, alpha=0.1) +
  scale_x_continuous(trans=log2_trans()) +
  labs(x="", y="") +
  facet_grid(.~ np.ratio) +
  theme(legend.position="none")

supp.all <- 
  ggplot(out.all, aes(x=lambda,y=error.sup, colour=method,group=method,lty=method)) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth", width=3, alpha=0.1) +
  scale_x_continuous(trans=log2_trans()) +
  labs(x="", y="") +
  facet_grid(.~ np.ratio) +
  theme(legend.position=c(.875, .74))

dprec <- qplot(factor(method), prec, data=out.all, group=method, colour=method, geom="boxplot", log="y") +
  labs(title=paste("Numerical accuracy, n=",n0), x="", y="Distance to optimum") +
  theme(legend.position="none")

dtime <- qplot(factor(method), times, data=out.all, group=method, colour=method, geom="boxplot", log="y") +
  labs(x="", y="Timing (sec.)")+
  theme(legend.position="none")

pdf(file="test,support.pdf", width=7, height=7)
##multiplot(test1, supp1,
##          test2, supp2,
##          test3, supp3, cols=2)
multiplot(test.all, supp.all, cols=1)
dev.off()

