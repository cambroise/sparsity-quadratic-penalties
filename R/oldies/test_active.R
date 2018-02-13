rm(list=ls())
dyn.load("functions/src/active_set_cpp.so")
active_set <- function(beta, xtx, A, x, D, xty, nabla.f, lambda, eps, solver="quadra") {

  solver_cpp <- switch(solver, 
                       quadra   = 0,
                       pathwise = 1,
                       fista    = 2)
  
  out <- .Call("active_set_cpp",
               beta    ,
               xtx     ,
               x       ,
               D       ,
               xty     ,
               nabla.f ,
               lambda  ,
               eps     ,
               solver_cpp)
  return(list(beta.A  = out$beta.A ,
              xtx.A   = out$xtx.A  ,
              nabla.f = out$nabla.f,
              active  = out$active ,
              gap     = out$gap    ,
              conv    = out$conv   ,
              iter    = c(out$iter)))
}

p <- 5
n <- 10
x <- matrix(rnorm(p*n),n,p)
y <- rnorm(10)

x.bar <- colMeans(x)
x     <- scale(x,x.bar,FALSE) 
y.bar <- mean(y)
y <- y-y.bar
norm <- sqrt(colSums(x^2))
x    <- scale(x, FALSE, norm)
xty <- c(crossprod(y,x))

gamma <- 0
D <- matrix(0,p,p)
diag(D) <- 1
D <- gamma * D

beta    <- numeric(0)
xtx     <- numeric(0)
active  <- numeric(0)
nabla.f <- -xty
lambda <- 1
eps <- 1e-5

out_quad <- active_set(beta, xtx, active, x, D, xty, nabla.f, lambda, eps, solver="quadra")
out_prox <- active_set(beta, xtx, active, x, D, xty, nabla.f, lambda, eps, solver="fista")
out_path <- active_set(beta, xtx, active, x, D, xty, nabla.f, lambda, eps, solver="pathwise")

source("functions/simulations.R")

n <- 3000
beta <- rep(c(0,2,0,-2),each=5)
beta <- rep(beta,100)
p <- length(beta)

cat("\nData generation")
data <- example_c_enet(n, beta, cor=0.25)
x <- data$x
y <- data$y

x.bar <- colMeans(x)
x     <- scale(x,x.bar,FALSE) 
y.bar <- mean(y)
y <- y-y.bar
norm <- sqrt(colSums(x^2))
x    <- scale(x, FALSE, norm)
xty <- c(crossprod(y,x))

gamma <- 0
D <- matrix(0,p,p)
diag(D) <- 1
D <- gamma * D

beta    <- numeric(0)
xtx     <- numeric(0)
active  <- numeric(0)
nabla.f <- -xty
lambda <- 0.5
eps <- 1e-5

time.quad <- system.time(out_quad <- active_set(beta, xtx, active, x, D, xty, nabla.f, lambda, eps, solver="quadra"))
time.prox <- system.time(out_prox <- active_set(beta, xtx, active, x, D, xty, nabla.f, lambda, eps, solver="fista"))
time.path <- system.time(out_path <- active_set(beta, xtx, active, x, D, xty, nabla.f, lambda, eps, solver="pathwise"))


