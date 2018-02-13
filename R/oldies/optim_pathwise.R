dyn.load("functions/src/pathwise.so")
pathwise <- function(x0, xtx, xty, lambda, eps) {
  return(.C("PATHWISE",
            as.integer(length(x0)),
            as.double(xtx),
            as.double(xty),
            as.double(lambda),
            as.double(diag(as.matrix(xtx))),
            as.integer(5000), # max # of iteration
            as.double(eps/length(x0)),
            x0 = as.double(x0),
            xk = as.double(x0),
            i  = as.integer(0))) 
}

dyn.load("functions/src/pathwise_cpp.so")
pathwise.cpp <- function(x0, xtx, xty, lambda, eps) {
  out <- .Call("pathwise_cpp", x0, as.matrix(xtx), xty, lambda, eps)
  return(list(xk=c(out$xk),i=out$i))
}
