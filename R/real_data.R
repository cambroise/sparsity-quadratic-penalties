rm(list=ls())
gc()
setwd("~/svn/sandbox/R/")
source("functions/simulations.R")
source("functions/plots.R")

## COLOR BLIND PALETTE
library(ggplot2)
library(scales)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scale_colour_manual(values=cbPalette)

## AVAILABLE DATASET
studies <- c("golub", "guedj", "lung", "julia", "prostate")
study <- "guedj"
data  <- getData(study)

lambda <- switch(study,
                 "golub"    = 10^seq(1, -1, length=50),
                 "guedj"    = 10^seq(1, -1, length=50),
                 "lung"     = 10^seq(1, -1, length=50),
                 "julia"    = 10^seq(1, -2, length=50),
                 "prostate" = 10^seq(1, -2, length=50))

## METHODS COMPARISON
thres.glmn <- 10^c(-12, -6, -1)
thres.spam <- 10^c(-6, -3, -1)
## main function
out <- average.test.error(K=100, data$x, data$y, lambda, thres.glmn, thres.spam, mycode=FALSE, gamma=0)

## Graphical outputs
out.fista <- out[out$methods != "glmn", ]
out.glmn  <- out[out$methods != "spam", ]

spam <-
  ggplot(out.fista, aes(x=lambda, y=error.mse, colour=methods, group=instance, lty=instance)) +
  geom_point(alpha=.3, size=1) +
  geom_smooth(alpha=.2, size=1) +
  scale_y_continuous(trans=log2_trans()) +
  scale_x_continuous(trans=log2_trans()) +
  opts(title="Estimated test error, quadrupen vs. SPAMs")

glmn <- 
ggplot(out.glmn, aes(x=lambda, y=error.mse, colour=methods, group=instance, lty=instance)) +
  geom_point(alpha=.3, size=1)+
  geom_smooth(alpha=.2, size=1) +
  scale_y_continuous(trans=log2_trans()) +
  scale_x_continuous(trans=log2_trans()) +
  opts(title="Estimated test error, quadrupen vs. glmnet")

multiplot(spam, glmn, cols=2)
