bfgs <- function(x0, xtx, xty, lambda, eps) {
  res <- optim(x0, method="BFGS",
               fn     = L,
               gr     = dL,
               xtx    = xtx,
               xty    = xty,
               lambda = lambda,
               control= list(reltol=(eps/length(x0))^2))
  i <- res$counts[1]
  names(i) <- NULL
  return(list(i=i,xk=res$par))
}

L  <- function(beta, xtx, xty, lambda) {
  .5*crossprod(beta,crossprod(xtx,beta))-crossprod(beta,xty)+lambda*sum(abs(beta))
}

dL <- function(beta, xtx, xty, lambda) {
  
  dL      <- rep.int(0,length(beta))
  nabla.f <- drop(crossprod(xtx,beta)) - xty
  n.nabla.f <- abs(nabla.f)

  ## add the elastic-net term in a second time...
  nabla.f <- nabla.f
  n.beta <- abs(beta)
  zero <- n.beta == 0
  
  ## terms whose norm1 is different from zero
  dL[!zero] <- nabla.f[!zero] + lambda * beta[!zero]/n.beta[!zero]
  
  ## terms whose norm2 is zero
  shrunken  <- zero & n.nabla.f > lambda
  
  ## shrink toward zero
  shrink <- lambda/n.nabla.f[shrunken]
  dL[shrunken] <- nabla.f[shrunken] * (1 - shrink)
  
  return(dL) 
}
