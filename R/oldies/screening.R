screening <- function(x,y, lambda.min=0.1, n.lambda=100, method="SAFE", penalty = "elastic.net", gamma=0) {

  screen.shot <- function(l) {
    apply( abs(xty) < l - sqrt(sum(y^2)) * sqrt(1 + gamma) * (lambda.max-l)/lambda.max, 1, which)
  }
  
  ## predictor are assumed to be normalized
  y <- y - mean(y)
  x <-  scale(x,TRUE,FALSE)
  ## unitary norm for each column of x
  norm <- sqrt(drop(colSums(x^2)))
  x    <- scale(x, FALSE, norm)
  xty <- crossprod(y,x)
  
  lambda.max <- max(abs(xty))
  lambda <- 10^seq(log10(lambda.max), log10(lambda.max*lambda.min) ,length=n.lambda)
  
  omitted <- sapply(lambda, screen.shot)

  return(list(omitted=omitted, lambda=lambda))
}
