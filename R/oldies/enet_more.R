rm(list=ls())
gc()
library(elasticnet)
library(glmnet)
library(R.utils)
sourceDirectory("functions")

seed <- sample(1:10000,1)
set.seed(seed)
## seed 6979 is nice

n <- 400
beta <- rep(c(0,2,0,-2),each=5)
beta <- rep(beta,10)
p <- length(beta)

cat("\nData generation")
data <- example_c_enet(n, beta, cor=0.25)
x <- data$x
y <- data$y

cat("\nGlmnet...")
time.glmnet <- system.time(out.glmnet <- glmnet(x, y, thres=1e-10))

gamma <- 0

cat("\nLarsen...")
## time.enet <- system.time(out.enet <- lars(x, y, use.Gram=FALSE))
time.enet <- system.time(out.enet <- lars(x, y))
##time.enet <- system.time(out.enet <- enet(x, y, lambda=gamma))

cat("\nBFGS...")
time.bfgs <- system.time(out.bfgs   <- crafter(x, y, gamma=gamma,   lambda =  sqrt(n) * out.glmnet$lambda,
                                                  optim.method="bfgs"))

cat("\nQuadra (full C++)...")
time.quad <- system.time(out.quad <- crafter.cpp(x, y,
                                                 gamma=gamma,
                                                 lambda =  sqrt(n) * out.glmnet$lambda))

coef.enet <- out.enet$beta[-length(out.enet$lambda),]
enet.fit <- new("fit",
                coefficients = Matrix(coef.enet),
                lambda       = out.enet$lambda,
                intercept    = rep(out.enet$mu,nrow(coef.enet)),
                monitoring   = list(),
                call         = out.enet$call)

coef.glmn <- t(out.glmnet$beta)
glmn.fit <- new("fit",
                coefficients = coef.glmn,
                lambda       = out.glmnet$lambda*sqrt(n),
                intercept    = out.glmnet$a0,
                monitoring   = list(),
                call         = out.glmnet$call)



d.quad <- dual_gap(out.quad,x,y)
d.bfgs <- dual_gap(out.bfgs,x,y)
d.enet <- dual_gap(enet.fit,x,y)
d.glmn <- dual_gap(glmn.fit,x,y)


par(mfrow=c(1,3))
plot(enet.fit, xvar="fraction", main= "enet")
plot(out.quad, xvar="fraction", main= "quadra")
boxplot(d.quad,
        d.bfgs,
        d.enet,
        d.glmn, ylab="dual gap",
        names=c("quadra","bfgs","larsen","glmnet"), las=3)

cat("\n TIMING\n")
time <- c(time.enet[3],time.bfgs[3],time.glmnet[3],time.quad[3])
names(time) <- c("Larsen","BFGS","Glmnet","Quadra") 
print(time)
