quadra_cpp <- function(x0, xtx, xty, lambda, eps) {

  A <- which(abs(x0) != 0)
  
  sgn    <- sign(xty)
  sgn[A] <- sign(x0[A])

  C <- chol.default(xtx)
  beta <- backsolve(C, backsolve(C, xty - lambda*sgn, transpose=TRUE))
  
  swap <- A[sign(beta[A]) != sgn[A]]

  ## if any sign swaps
  if (any(swap)) {
    swapped <- beta[swap]
    ## first, go to zero for the swapped variables 
    gam <- -x0[swap]/(beta-x0)[swap]
    beta <- x0 + min(gam)*(beta-x0)
    ## second, solve the problem after swaping the signs of the
    ## increminated variables
    x0 <- beta
    x0[swap] <- -swapped
    A <- which(abs(x0) != 0)
    sgn    <- sign(xty)
    sgn[A] <- sign(x0[A])
    beta.tilde <- backsolve(C, backsolve(C, xty - lambda*sgn, transpose=TRUE))
  } else {
    beta.tilde <- NULL
  }
 
  return(list(xtilde = beta.tilde, xk = beta, i=ifelse(any(swap),2,1), swap = swap))
}

dyn.load("functions/src/quadra_cpp.so")
quadra.cpp <- function(x0, xtx, xty, lambda, eps) {
  out <- .Call("quadra_cpp", x0, as.matrix(xtx), xty, lambda, eps)
  if (is.null(out$swap)) {
    return(list(xtilde=NULL, xk=c(out$xk),i=1, swap= integer(0)))
  } else {
    return(list(xk=c(out$xk),i=2,xtilde=c(out$xtilde),swap=c(out$swap+1)))
  }
}
