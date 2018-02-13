rm(list=ls())
library(Matrix)
dyn.load("sp_prod.so")

n <- 1000
p <- 5000
prob <- 0.01
size <- 5

Xd  <- matrix(rbinom(n*p, size, prob),n,p)
X   <- as(Xd, "dgTMatrix")
SX  <- list(Xi = X@i, Xj = X@j, Xx = X@x)
X2  <- as(Xd, "dgCMatrix")
SX2 <- list(Xi = X2@i, Xp = X2@p, Xnp = diff(X2@p), Xx = X2@x)

j <- which(colSums(X2) !=0)[1]

t.Matrix <- system.time(out.Matrix <- t(X) %*% X[, j])
t.matrix <- system.time(out.matrix <- t(Xd) %*% Xd[, j])
t.mine   <- system.time(out.mine   <- .Call("sp_prod" , SX , j-1, n, p))
t.mine2  <- system.time(out.mine2  <- .Call("sp_prod2", SX2, j-1, nrow(X2), ncol(X2)))

