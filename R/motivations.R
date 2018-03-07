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
set.seed(seed)
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
lambda <- 10^seq(0.5, -1, len=100)

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
}

library(dplyr)

dplot.all <- out.all %>% 
  mutate(rmse = sqrt(error.mse), 
          precision = ifelse(tp == 0, 0, tp /(tp + fp)),
          recall = ifelse(tp == 0, 0, tp /(tp + fn))) %>% 
  mutate(f_measure = ifelse ( (precision == 0 | recall == 0), 0,  2 * (precision * recall) / (precision + recall))) %>% 
  rename(hamming = error.sup, timings = times, accuracy = prec) %>% 
  select(-error.cla, -error.mse, -steps, -package) %>% 
  filter(method != "glmnet high")

dplot <- dplot.all %>% select("lambda", "method", "np.ratio", "precision", "recall", "rmse") %>% gather(key = "measure", value = "value", "rmse", "precision", "recall")

p <- ggplot(dplot, aes(x = lambda, y = value, colour = method, group = method, lty = method)) +
  stat_summary(fun.data="mean_cl_normal", geom = "smooth", alpha = 0.1) +
  coord_trans(x = "log10") + labs(x="", y="") +
  facet_grid(measure ~ np.ratio, scales = "free") + theme(legend.position="none") + theme_minimal()

ggsave(filename = "../figures/accuracy.pdf", width = 10, height = 7)
