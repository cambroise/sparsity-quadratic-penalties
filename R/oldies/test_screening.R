rm(list=ls())
setwd("~/svn/sandbox/R")
library(scoop)
library(R.utils   , quietly = TRUE)
gc()
sourceDirectory("functions")

n <- 200
p <- 2000
grp <- rep(1:100,each=20)
gamma <- c(0.001,0.01,0.1,1,10)

data <- opti_book_settings(n, p, 0.1, "low")

screen <- screening(data$x,data$y,gamma=1,lambda.min=0.1)
plot(screen$lambda,sapply(screen$omitted,length))
