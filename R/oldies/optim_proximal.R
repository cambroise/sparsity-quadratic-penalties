dyn.load("functions/src/fista.so")
fista <- function(x0, xtx, xty, lambda, eps) {
  return(.C("FISTA",
            as.integer(length(x0)),
            as.double(xtx),
            as.double(xty),
            as.double(lambda),
            as.double(max(eigen(xtx,TRUE,TRUE)$values)),
            as.integer(5000), # max # of iteration
            as.double(eps/length(x0)),
            x0 = as.double(x0),
            xk = as.double(x0),
            i  = as.integer(0)))
}

dyn.load("functions/src/fista_cpp.so")
fista.cpp <- function(x0, xtx, xty, lambda, eps) {
  out <- .Call("fista_cpp", x0, as.matrix(xtx), xty, lambda, eps)
  return(list(xk=c(out$xk),i=out$i))
}

