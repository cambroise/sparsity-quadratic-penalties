crafter <- function(x,
                    y,
                    gamma        = 10,
                    D            = NULL,
                    intercept    = TRUE,
                    normalize    = TRUE,
                    lambda       = NULL,
                    n.lambda     = ifelse(is.null(lambda),100,length(lambda)),
                    lambda.min   = ifelse(gamma == 0,1e-1,1e-3),
                    verbose      = FALSE,
                    eps          = 1e-3,
                    max.features = min(ncol(x), 1000),
                    max.iter     = 2*max.features,
                    optim.method = "quadra",
                    profiling    = TRUE) {
  
  solver <- switch(optim.method,
                   bfgs         = bfgs,
                   fista        = fista,
                   quadra       = quadra,
                   pathwise     = pathwise,
                   fista.cpp    = fista.cpp,
                   quadra.cpp   = quadra.cpp,
                   pathwise.cpp = pathwise.cpp)

  if(profiling)
    prof.out <- paste(c("profiling_"),optim.method,".txt",sep="")
  
  ## ======================================================
  ## PRETREATMENT
  n <- nrow(x)
  p <- ncol(x)
  if (is.null(colnames(x))) {colnames(x) <- 1:p}
  ## centering the data
  if (intercept) {
    x.bar <- colMeans(x)
    x     <- scale(x,x.bar,FALSE) 
    y.bar <- mean(y)
    y <- y-y.bar
  }
  ## normalizing the data
  if (normalize) {
    norm <- sqrt(colSums(x^2))
    x    <- scale(x, FALSE, norm)
  } else {
    norm <- rep(1,p)
  }
  ## defaut structure is... none (identity matrix)
  if (is.null(D)) {
    D <- matrix(0,p,p)
    diag(D) <- 1
  }
  D <- gamma * D
  
  ## ========================================================
  ## INITIALIZATION OF THE MAIN OBJECTS 
  Beta    <- c() ## parameters
  beta.A  <- c() ## active parameters
  xty     <- c(crossprod(y,x))
  xtx.A   <- c()
  mu      <- c()
  nabla.f <- -xty ## smooth part of the gradient
  active  <- c()  ## the successively activated variables
  monitoring  <- list(it.active=c(), it.optim=c(), dual.gap=c(), timing=c())
  ends <- c("converged", "max # of iterates reached", "no improvement...")
  
  ## ========================================================
  ## LOOP OVER LAMBDA 
  ## Generate a grid of lambda if none has benne provided
  if (is.null(lambda)) {
    lmax <- max(abs(xty))
    lambda <- 10^seq(from=log10(lmax),to=log10(lambda.min), length=n.lambda)
  }
  m <- 0
  ## Starting clock
  if (profiling) {Rprof("tmp.out")}
  t0 <- proc.time()
  for (l in lambda) {
    if (verbose) cat("\nlambda =",l)
    m    <- m + 1
    iter <- c()
    time <- c()
    ## _____________________________________________________________
    ##
    ## START THE ACTIVE SET ALGORITHM
    ## _____________________________________________________________
    ##
    dual.all <- pmax.int(0,abs(nabla.f) - l)
    if (any(active))
      dual.all[active] <- abs(nabla.f[active] + l * sign(beta.A))
    gap  <- max(dual.all)
    conv <- c(gap < eps, FALSE, FALSE)
    while (!any(conv)) {
      ## _____________________________________________________________
      ##
      ## (1) VARIABLE ACTIVATION IF APPLICABLE
      ## _____________________________________________________________

      new <- which.max(dual.all)
      if (!(new %in% active)) {
        if (verbose) cat(" +",new)
        active <- c(active, new)
        xtx.A  <- cbind(xtx.A,crossprod(x,x[,new]) + D[,new])
        beta.A <- c(beta.A,0)
      }
      ## _____________________________________________________________
      ##
      ## (2) OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
      ## _____________________________________________________________
      ##

      res <- solver(beta.A, xtx.A[active,], xty[active], l, eps)
      iter <- c(iter,res$i)
      ## _____________________________________________________________
      ##
      ## (3) VARIABLE DELETION IF APPLICABLE
      ## _____________________________________________________________
      ##

      if (is.null(res$xtilde)) {
        beta.A <- res$xk
        nabla.f      <- -xty + tcrossprod(beta.A, xtx.A)
        ## unactivate variables with null gradient
        null <- which(abs(beta.A) < eps & abs(nabla.f[active]) < l + eps)
        if (length(null) > 0) {
          if (verbose) cat(" - ",active[null])
          xtx.A  <- xtx.A[,-null]
          active <- active[-null]
          beta.A <- beta.A[-null]
        }
      } else {
        ## variables with sign swaps (only relevent for the quadra solver)
        ## two solutions: one with zeroed variables after sign swap,
        ## and the other one with sign effectively swapped
        ## We keep the one with the smallest dual gap.
        nabla.f1 <- -xty + tcrossprod(res$xk    , xtx.A)
        nabla.f2 <- -xty + tcrossprod(res$xtilde, xtx.A)
        dual1 <- abs(nabla.f1[active[res$swap]]) - l
        dual2 <- abs(nabla.f2[active[res$swap]]  + l*sign(res$xtilde[res$swap]))

        if (min(dual2) < eps) {
          beta.A <- res$xtilde
          nabla.f      <- nabla.f2
        } else {
          ## in this case, the variable is unactivated
          beta.A <- res$xk
          if (verbose) cat(" - ",active[res$swap])
          xtx.A  <- xtx.A[,-res$swap]
          beta.A <- beta.A[-res$swap]
          active <- active[-res$swap]
          nabla.f      <- nabla.f1
        }
      }

      ## _____________________________________________________________
      ##
      ## (4) OPTIMALITY TESTING
      ## _____________________________________________________________
      ##
      dual.all <- pmax.int(0,abs(nabla.f) - l)
      if (any(active))
        dual.all[active] <- abs(nabla.f[active] + l * sign(beta.A))
      gap.old <- gap
      gap  <- max(dual.all)
      conv <- c(gap < eps, length(iter) == max.iter, abs(gap-gap.old) < eps)
      
    } ## END THE ACTIVE SET HERE
    
    ## Stop if the maximal number of features requested has been reached
    if (length(active) > max.features) {
      break
    }
    
    ## When the algorithm does not converge, send a message
    if (which(conv)[1] != 1) {
      if (verbose) {
        cat("\n",ends[conv], "with max.min||subgradient|| =", gap)
      }
      time <- NA
    } else {
      time <- (proc.time()-t0)[3]
    }
    names(time) <- round(l,2)
    monitoring$it.active <- c(monitoring$it.active, length(iter))
    monitoring$it.optim  <- c(monitoring$it.optim , iter)
    monitoring$dual.gap  <- c(monitoring$dual.gap , gap)
    monitoring$timing    <- c(monitoring$timing   , time)

    ## Move to the next lambda
    beta <- rep(0,p)
    if (normalize) {
      beta[active] <- beta.A/norm[active]
    } else {
      beta[active] <- beta.A
    }
    Beta <- rbind(Beta,beta)
    ## intercept
    if (intercept) {
      mu   <- c(mu,y.bar - crossprod(beta.A/norm[active]),x.bar[active])
    }
    
    ## Stop if the path of solutions has stabilized
    ##     if (length(active) > 0 & m > 1) {
    ##       if (sqrt(sum((Beta[(m-1),]-Beta[m, ])^2)) < eps) {
    ##         break
    ##       }
    ##     }
    
  } ## END OF THE LOOP OVER LAMBDA
  if (profiling) {
    Rprof(NULL)
    system(paste(c("R CMD Rprof tmp.out > "),prof.out))
    system("rm *.out")
  }
  
  ## fit's construction
  return(new("fit",
             coefficients = Matrix(Beta),
             intercept    = ifelse(is.null(mu),0,mu),
             lambda       = lambda[1:m],
             monitoring   = monitoring,
             call         = match.call()))

}
