library(quadrupen)
source("functions/simulations.R")

seed <- sample(1:10000,1)
set.seed(seed)

sample.scenario <- "low"
sparse.scenario <- "high"
eps <- 1e-2

## Sparsity settings
beta <- switch(sparse.scenario,
               "low"  = rep(c(2,-2,0),c(20,20,60)),
               "med"  = rep(c(2,-2,0),c(15,15,70)),
               "high" = rep(c(2,-2,0),c(5,5,90)))
p <- length(beta)

## Dimensional settings for Lasso run
n <- switch(sample.scenario, "lower" = p/4, "low" = p/2, "med"= p, "high" = 2*p)
l <- switch(sample.scenario, "lower" = 0.001, "low" = 0.001, "med"= 0.0005, "high" = 0.0001)
file.out <- paste("output/compa-optim_p=",p,",s=",sum(beta!=0),",n=",n,",eps=",eps,".RData", sep="")

nsim <- 5
nlambda <- 50
ngamma  <- 50
gammas <- 10^seq(-2,1,len=ngamma)

## SMALL CORRELATION
cat("\n\n Generating small correlations data set...")
simu.small <- sim(nsim,n,beta,cor=0.1)
cat("\n\n Generating medium correlations data set...")
simu.medi <- sim(nsim,n,beta,cor=0.5)
cat("\n\n Generating large correlations data set...")
simu.larg <- sim(nsim,n,beta,cor=0.8)

lambda.max <- max(simu.small$lambda.max,simu.medi$lambda.max,simu.larg$lambda.max)
lambdas <- 10^seq(from=log10(lambda.max),to=log10(lambda.max*l), length=nlambda)

cat("\n\n Computation for small correlations...")
cat("\n Crafter...")
res.quad.small <- run(simu.small, lambdas, gammas, method="crafter",eps=eps, mc.cores=mc.cores)
cat("\n Pathwise...")
res.path.small <- run(simu.small, lambdas, gammas, method="pathwise",eps=eps, mc.cores=mc.cores)
cat("\n Proximal...")
res.prox.small <- run(simu.small, lambdas, gammas, method="fista",eps=eps, mc.cores=mc.cores)

cat("\n\n Computation for medium correlations...")
cat("\n Crafter...")
res.quad.medium <- run(simu.medi, lambdas, gammas, method="crafter",eps=eps, mc.cores=mc.cores)
cat("\n Pathwise...")
res.path.medium <- run(simu.medi, lambdas, gammas, method="pathwise",eps=eps, mc.cores=mc.cores)
cat("\n Proximal...")
res.prox.medium <- run(simu.medi, lambdas, gammas, method="fista",eps=eps, mc.cores=mc.cores)

cat("\n\n Computation for large correlations...")
cat("\n Crafter...")
res.quad.large <- run(simu.larg, lambdas, gammas, method="crafter",eps=eps, mc.cores=mc.cores)
cat("\n Pathwise...")
res.path.large <- run(simu.larg, lambdas, gammas, method="pathwise",eps=eps, mc.cores=mc.cores)
cat("\n Proximal...")
res.prox.large <- run(simu.larg, lambdas, gammas, method="fista",eps=eps, mc.cores=mc.cores)

rm(simu.larg, simu.medi, simu.small)
save.image(file=file.out)
