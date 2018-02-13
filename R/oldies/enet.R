rm(list=ls())
gc()
library(elasticnet)
library(R.utils)
sourceDirectory("functions")

seed <- sample(1:10000,1)
set.seed(seed)

n <- 100
beta <- rep(c(0,2,0,-2),each=5)
beta <- rep(beta,10)
p <- length(beta)

data <- example_c_enet(n, beta, cor=0.25)
x <- data$x
y <- data$y

gamma <- 0

time.enet <- system.time(out.enet <- enet(x, y, lambda=gamma))

cat("\nBFGS...")
time.bfgs <- system.time(out.bfgs   <- crafter(x, y, gamma=gamma,
                                                 optim.method="bfgs"))
cat("\nFISTA (C++)...")
time.prox <- system.time(out.prox  <- crafter.cpp(x, y, gamma=gamma, 
                                                  optim.method="fista"))

cat("\nPATHWISE (C++)...")
time.path <- system.time(out.path  <- crafter.cpp(x, y, gamma=gamma, 
                                                 optim.method="pathwise"))

cat("\nQuadra (full C++)...")
time.quad <- system.time(out.quad <- crafter.cpp(x, y, gamma=gamma,
                                                 optim.method="quadra"))
repo <- paste("runtimes_gamma=",gamma,"_np=",n,"x",p,sep="")
system(paste("mkdir",repo))
system(paste("mv profiling*",repo))

coef.enet <- out.enet$beta[-length(out.enet$penalty),]/(1+gamma)
enet.fit <- new("fit",
                coefficients = Matrix(coef.enet),
                lambda       = out.enet$penalty[-length(out.enet$penalty)]/2,
                intercept    = rep(out.enet$mu,nrow(coef.enet)),
                monitoring   = list(),
                call         = out.enet$call)

par(mfrow=c(1,3))
plot(enet.fit, xvar="fraction", main= "enet")
plot(out.quad, xvar="fraction", main= "quadra")
boxplot(out.bfgs@monitoring$dual.gap,
        out.prox@monitoring$dual.gap,
        out.path@monitoring$dual.gap,
        out.quad@monitoring$dual.gap, ylab="dual gap",
        names=c("bfgs","proximal","pathwise","quadra"), las=3)
d <- sum( (out.quad@coefficients - out.bfgs@coefficients )^2 )
cat("\n ||quadra-bfgs||^2 =",d,"\n")

time <- c(time.enet[3],time.bfgs[3],time.prox[3],time.path[3],time.quad[3])
names(time) <- c("Larsen","BFGS","FISTA","Pathwise","Quadra") 
print(time)
