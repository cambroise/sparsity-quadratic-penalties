library(parallel)

rmvnorm.const.cor <- function(n, p, rho=0.5) {

  C <- chol.const.cor(p,rho)
  rmv <- function() {
    Z <- rnorm(p)
    return(C$cii*Z + c(0,cumsum(Z[-p]*C$c.i)))
  }
  
  return(t(replicate(n, rmv())))
}

chol.const.cor <- function(p,rho) {

  q <- 0     ## sum_k^(i-1) c.k^2
  cii <- 1   ## current value of the diagonal 
  c.i <- rho ## current value of the off-diagonal term
  Cii <- cii ## vector iof diagonal terms
  C.i <- c.i ## vector of off-diagonal terms
  for (i in 2:p) {
    q <- q + c.i^2
    cii <- sqrt(1-q)
    c.i <- (rho-q)/cii
    Cii <- c(Cii,cii)
    C.i <- c(C.i,c.i)
  }
  return(list(cii=Cii, c.i=C.i[-p]))
}

rlm <- function(n, beta, corr=0.5, mu=3, r2=0.8, sparse=FALSE,prob=0.01,size=5) {
  p <- length(beta)
  Sigma  <- matrix(corr, p,p)
  diag(Sigma) <- 1
  if (sparse) {
    X <- matrix(rbinom(n*p, size, prob),n,p)
  } else {
    X <- rmvnorm.const.cor(n,p,corr)
  }
  sigma  <- sqrt((1-r2)/r2 * t(beta) %*% Sigma %*% beta)
  epsilon <- rnorm(n) * sigma
  y <- mu + X %*% beta + epsilon
  r2 <- 1 - sum(epsilon^2) / sum((y-mean(y))^2)
  
  return(list(y=y, x=X, r2=r2))
}

rlm.IC.lasso <- function(n,beta,sigma,mu=3,n.wish=p) {
  require(bayesm)
  require(mvtnorm)
  p <- length(beta)
  S <- beta != 0
  IC <- FALSE

  while (!IC) {
    W <- cov2cor(rwishart(n.wish,diag(p))$W)
    X <- rmvnorm(n, sigma=W)
    Gram <- t(X) %*% X/n
    eta <- max(Gram[!S,S] %*% solve(Gram[S,S]) %*% sign(beta[S]))
    IC <- eta >= 1
  }
  epsilon <- rnorm(n) * sigma
  y <- mu + X %*% beta + epsilon
  r2 <- 1 - sum(epsilon^2) / sum((y-mean(y))^2)
  return(list(y=y, x=X, r2=r2, IC = IC, eta=eta))
}

run <- function(simu, lambdas, gammas, method="crafter", eps= 1e-3, verbose=FALSE, mc.cores=4) {
  
  nsim   <- ncol(simu$data)
  fun <- function(gam,x,y) {
    craftout <- elastic.net(x,y,lambda2=gam,lambda1=lambdas, 
                            control = list(method=method, thres=eps, verbose=verbose))
    craftout@monitoring$pensteps.timer[craftout@monitoring$converge != "converged"] <- NA
    craftout@monitoring$it.optim[craftout@monitoring$converge != "converged"] <- NA

    timer <- craftout@monitoring$pensteps.timer
    qiopt <- quantile(craftout@monitoring$it.optim)
    qiact <- quantile(craftout@monitoring$it.active)
    
    return(list(timer=timer,qiopt=qiopt, qiact=qiact))
  }
  
  res.timer <- matrix(0, length(lambdas),length(gammas))
  nor.timer <- matrix(0, length(lambdas),length(gammas))
  res.qiopt <- matrix(0,5,length(gammas))
  nor.qiopt <- matrix(0,5,length(gammas))
  res.qiact <- matrix(0,5,length(gammas))
  nor.qiact <- matrix(0,5,length(gammas))

  colnames(res.timer) <- round(gammas , 2)
  rownames(res.timer) <- round(lambdas, 2)
  colnames(res.qiopt) <- round(gammas , 2)
  colnames(res.qiact) <- round(gammas , 2)
  
  cat("\n===========================================")
  cat("\n# simu  (over ",nsim,") : ", sep="")
  
  for (s in 1:nsim) {
    cat("",s)

    out   <- do.call(cbind, mclapply(gammas, fun, x=simu$data[,s]$x, y=simu$data[,s]$y, mc.cores=mc.cores))
    timer <- matrix(unlist(out[1, ]), ncol=length(gammas))
    qiopt <- matrix(unlist(out[2, ]), ncol=5)
    qiact <- matrix(unlist(out[3, ]), ncol=5)
    
    res.timer[!is.na(timer)] <- res.timer[!is.na(timer)] + timer[!is.na(timer)]
    nor.timer[!is.na(timer)] <- nor.timer[!is.na(timer)] + 1
    
    res.qiopt[!is.na(qiopt)] <- res.qiopt[!is.na(qiopt)] + qiopt[!is.na(qiopt)]
    nor.qiopt[!is.na(qiopt)] <- nor.qiopt[!is.na(qiopt)] + 1

    res.qiact[!is.na(qiact)] <- res.qiact[!is.na(qiact)] + qiact[!is.na(qiact)]
    nor.qiact[!is.na(qiact)] <- nor.qiact[!is.na(qiact)] + 1
  }
  return(list(timer=res.timer/nor.timer, qiopt=res.qiopt/nor.qiopt, qiact=res.qiact/nor.qiact)) 
  
}

sim <- function(nsim, n, beta, cor, nlambda, lambda.min) {
  
  ## DATA GENERATION

  data <- replicate(nsim, rlm(n, beta, cor))
  out <- sapply(1:nsim, function(s) {
    ## scenario
    x <- data[,s]$x
    y <- data[,s]$y
    ## normalizing the data
    x     <- scale(x,colMeans(x),FALSE) 
    norm <- sqrt(colSums(x^2))
    x    <- scale(x, FALSE, norm)
    y <- y-mean(y)
    return(c(max(abs(crossprod(y,x))), data[,s]$r2))
  })
    
  return(list(data=data,lambda.max=max(out[1,]),r2.mean=mean(out[2,])))
}

glm2fit <- function(glm.obj,n) {
  coef.glmn <- t(glm.obj$beta)
  p <- ncol(coef.glmn)

    return(new("quadrupen",
                 coefficients  = Matrix(coef.glmn),
                 active.set    = Matrix(),
                 intercept     = TRUE       ,
                 mu            = glm.obj$a0  ,
                 normx         = rep(1,p)   ,
                 fitted        = Matrix()   ,
                 residuals     = Matrix()   ,
                 df            = 1        ,
                 r.squared     = 1        ,
                 penscale      = rep(1,p)   ,
                 penalty       = "lasso",
                 naive         = NULL    ,
                 lambda1       = glm.obj$lambda*sqrt(n),
                 lambda2       = 0  ,
                 control       = list()     ,
                 monitoring    = list()
                   ))
}

lar2fit <- function(lar.obj, xbar) {
  ## removing last lambda value, artificially added
  ## at the end of the LARS algorithm (saturated model)
  coef.lar <- lar.obj$beta[-length(lar.obj$lambda),]
  p <- ncol(coef.lar)
  intercept <- rep(lar.obj$mu,nrow(coef.lar))
  intercept <- c(intercept - coef.lar %*% xbar)
  lar.fit <- new("quadrupen",
                 coefficients  = Matrix(coef.lar),
                 active.set    = Matrix(),
                 intercept     = TRUE       ,
                 mu            = intercept  ,
                 normx         = rep(1,p)   ,
                 fitted        = Matrix()   ,
                 residuals     = Matrix()   ,
                 df            = 1        ,
                 r.squared     = 1        ,
                 penscale      = rep(1,p)   ,
                 penalty       = "lasso",
                 naive         = NULL    ,
                 lambda1       = lar.obj$lambda,
                 lambda2       = 0  ,
                 control       = list()     ,
                 monitoring    = list()
                 )
}

getlmax <- function(data) {
  x <- data$x
  y <- data$y
  x     <- scale(x,colMeans(x),FALSE) 
  norm <- sqrt(colSums(x^2))
  x    <- scale(x, FALSE, norm)
  y <- y-mean(y)
  return(max(abs(crossprod(y,x))))
}

lasso.spam <- function(x,y, lambda, type="lars", tol=1e-6, lambda2=0, lambda3=0) {

  if (type == "lars") {
    r.start <- proc.time()
    p <- ncol(x)
    n <- nrow(x)
    msteps <- min(n,p)
    ybar <- mean(y)
    xbar <- colMeans(x)
    normx <- sqrt(drop(colSums(scale(x,xbar,FALSE)^2)))
    mu <- ybar
    coef.spm <- t(spams.lasso(scale(y,ybar,FALSE), scale(x,xbar,normx),
                              return_reg_path=TRUE, lambda1=min(lambda),
                              max_length_path= msteps)[[2]])
    coef.spm  <- scale(coef.spm,FALSE,normx)
    intercept <- rep(mu,nrow(coef.spm))
    intercept <- c(intercept - coef.spm %*% xbar)
    timer <- (proc.time() - r.start)[3]
  }

  if (type == "fista") {
    r.start <- proc.time()

    p <- ncol(x)
    n <- nrow(x)
    msteps <- min(n,p)
    ybar <- mean(y)
    xbar <- colMeans(x)
    normx <- sqrt(drop(colSums(scale(x,xbar,FALSE)^2)))
    mu <- ybar
    coef.spm <- c()
    beta0 <- matrix(rep(0,p),p,1)
    for (l in lambda) {
      beta0 <- spams.fistaFlat(scale(y,ybar,FALSE), scale(x,xbar,normx), beta0,
                               loss="square", lambda1=l, lambda2=lambda2, regul="elastic-net", tol=tol)
      coef.spm <- cbind(coef.spm,beta0)
    }
    coef.spm  <- scale(t(coef.spm),FALSE,normx)
    intercept <- rep(mu,nrow(coef.spm))
    intercept <- c(intercept - coef.spm %*% xbar)
    
    timer <- (proc.time() - r.start)[3]
  }

    return(new("quadrupen",
                 coefficients  = Matrix(coef.spm),
                 active.set    = Matrix(),
                 intercept     = TRUE       ,
                 mu            = intercept  ,
                 normx         = rep(1,p)   ,
                 fitted        = Matrix()   ,
                 residuals     = Matrix()   ,
                 df            = 1        ,
                 r.squared     = 1        ,
                 penscale      = rep(1,p)   ,
                 penalty       = "lasso",
                 naive         = NULL    ,
                 lambda1       = lambda,
                 lambda2       = 0  ,
                 control       = list()     ,
                 monitoring    = list(external.timer=timer)
                   ))
  

}

test.lasso <- function(data, thres, fit.ref, package=c("glmnet","quadrupen","spams.lars","spams.prox")) {

  objc.ref <- objective(fit.ref,data$x,data$y)
  lambda   <- fit.ref@lambda1
    
  if ("glmnet" %in% package) {
    gc()
    time <- system.time(fit <- glmnet(data$x,data$y,
                                      thres=thres,
                                      lambda=lambda/sqrt(nrow(data$x)),
                                      pmax=min(ncol(data$x),nrow(data$x))))[3]

    objc <- objective(glm2fit(fit,nrow(data$x)),data$x,data$y)
    prec <- sqrt(sum( (objc - objc.ref)[-length(lambda)]^2)/length(lambda))

    if (fit$npasses > 1e5) {
      time <- NA
      prec <- NA
    }
  }
  if ("quadrupen" %in% package) {
    gc()
    fit  <- lasso(data$x,data$y, lambda, control=list(timer=TRUE))
    time <- fit@monitoring$external.timer
    objc <- objective(fit,data$x,data$y)
    prec <- sqrt(sum( (objc - objc.ref)[-length(lambda)]^2))/length(lambda)

    if (any(fit@monitoring$status != "converged")) {
      time <- NA
      prec <- NA
    }
  }

  if ("quadrupen.prox" %in% package) {
    gc()
    fit  <- lasso(data$x,data$y, lambda, control=list(timer=TRUE, method="fista", threshold=thres, max.iter=1000))
    time <- fit@monitoring$external.timer
    objc <- objective(fit,data$x,data$y)
    prec <- sqrt(sum( (objc - objc.ref)[-length(lambda)]^2))/length(lambda)

    if (any(fit@monitoring$status != "converged")) {
      time <- NA
      prec <- NA
    }
  }

  if ("quadrupen.path" %in% package) {
    gc()
    fit  <- lasso(data$x,data$y, lambda, control=list(timer=TRUE, method="pathwise", threshold=thres, max.iter=1000))
    time <- fit@monitoring$external.timer
    objc <- objective(fit,data$x,data$y)
    prec <- sqrt(sum( (objc - objc.ref)[-length(lambda)]^2))/length(lambda)

    if (any(fit@monitoring$status != "converged")) {
      time <- NA
      prec <- NA
    }
  }
  
  if ("spams.lars" %in% package) {
    gc()
    fit  <- lasso.spam(data$x,data$y,lambda,type="lars")
    time <- fit@monitoring$external.timer
    gc()
    objc <- objective(fit,data$x,data$y)
    prec <- sqrt(sum( (objc - objc.ref)[-length(lambda)]^2))/length(lambda)
  }

  if ("spams.prox" %in% package) {
    gc()
    fit  <- lasso.spam(data$x,data$y,lambda,type="fista",tol=thres)
    time <- fit@monitoring$external.timer
    objc <- objective(fit,data$x,data$y)
    prec <- sqrt(sum( (objc - objc.ref)[-length(lambda)]^2))/length(lambda)
  }
  
  res <- c(time, prec, objc)
  names(res) <- c("CPU time in sec.","Precision", "Objective")
  return(res)
}

test.enet <- function(gamma, data) {

  time.larsen  <- system.time(fit.lars <- enet(data$x,data$y,lambda=gamma))[3]

  time.crafter <- system.time({nlbd <- length(fit.lars$penalty)
                               lmax <- fit.lars$penalty[1]
                               lmin <- fit.lars$penalty[nlbd-1]
                               lambda <- 10^seq(from=log10(lmax), log10(lmin), len=nlbd)
                               fit.craf <- elastic.net(data$x,data$y,gamma=gamma,
                                                       lambda=lambda)})[3]
  
  return(c(time.larsen,time.crafter))
}

objective <- function(fit, x, y, gamma=0, mu=NULL) {
  
  J <- function(lambda) {
    l <- match(lambda, fit@lambda1)
    w <- fit@coefficients[l, ]
    if (is.null(mu)) {
        mu <- fit@mu[l]
    }
    1/length(y) * sum((y - x %*% w - mu)^2) + l * sum(abs(w)) + gamma * sqrt( sum((w)^2) )

  }
  
  return(sapply(fit@lambda1, J))
}

error <- function(fit, x.test, y.test, type="mse", thres=0.5) {

  y.hat <- predict(fit, x.test)

  if (type == "class") {
    y.hat[y.hat < thres] <- 0
    y.hat[y.hat >= thres] <- 1
    err <- apply(y.hat, 2, function(y) sum(y != y.test))
  }

  if (type == "mse") {
    err <- apply(y.hat, 2, function(y) mean((y - y.test)^2))
  }
  
  return(err)
}

outsample.test.error <- function(K=100, n, beta, cor, r2, lambda, thres.glmn, prop.test=1/4, mycode=FALSE, mc.cores=4) {

  cat("\n\nTEST SET SAMPLE NUMBER")

  one.simu <- function(i) {
      cat(" ",i)
      beta <- sample(beta)
      data <- rlm(n,beta,cor=cor,r2=r2)
      x <- data$x
      y <- data$y

      test <- sample(1:length(y),floor(length(y)*prop.test))
      train <- setdiff(1:length(y),test)
      return(compute_error(x[train, ], x[test, ], y[train], y[test], thres.glmn, thres.spam, beta, lambda, mine=mycode, gamma=gamma))
  }
  out <- do.call(rbind, mclapply(1:K, one.simu, mc.cores=mc.cores))
  sample <- factor(rep(1:K, each=nrow(out) %/% K))
  return(cbind(out, sample))

}

compute_error <- function(x.train, x.test, y.train, y.test, thres.glmn, thres.spam, true.beta, lambda,
                          max.feat=length(y.train)-1, thres.class=0.5, mine=FALSE, gamma=0) {
  require(glmnet)
##  require(spams)
  require(quadrupen)
  ## sample size
  n <- length(y.train)
  n.lambda <- length(lambda)
  
  ## quadratic solver
##  cat("\nQuadratic solver...")
  lasso.quadra <- lasso(x.train, y.train, lambda1=lambda, control=list(method="quadra", timer=TRUE))
  ##  cat(' ', length(lasso.quadra@lambda))

  ## computing error and timing
  error.quad.mse <- error(lasso.quadra, x.test, y.test)
  error.quad.cla <- error(lasso.quadra, x.test, y.test, type="class", thres=thres.class)
  error.quad.sup <- apply(lasso.quadra@coefficients, 1, function(x) sum(sign(x) != sign(true.beta)) )
  
  tp.quad <- apply(lasso.quadra@coefficients, 1, function(x) length(intersect(which(x != 0), which(true.beta != 0))))
  tn.quad <- apply(lasso.quadra@coefficients, 1, function(x) length(intersect(which(x == 0), which(true.beta == 0))))
  fp.quad <- apply(lasso.quadra@coefficients, 1, function(x) length(intersect(which(x != 0), which(true.beta == 0))))
  fn.quad <- apply(lasso.quadra@coefficients, 1, function(x) length(intersect(which(x == 0), which(true.beta != 0))))

  times.quad     <- rep(lasso.quadra@monitoring$external.timer, n.lambda)

  objc.ref <- objc <- objective(lasso.quadra,x.train,y.train)
  prec <- sqrt(sum(pmax(0,lasso.quadra@monitoring$max.grad)^2))/length(lambda)
  prec.quad <- rep(prec, n.lambda)
#  cat(" max.grad =", max(lasso.quadra@monitoring$max.grad))
#  cat(" min.grad =", min(lasso.quadra@monitoring$max.grad))
#  cat(" prec =",prec)
#  cat(" not converged =", sum(lasso.quadra@monitoring$status != "converged"))

##  cat("\nGlmnet for various thresholds...")
  ## glmnet for various thresholds...
  error.glmn.mse <- c()
  error.glmn.cla <- c()
  error.glmn.sup <- c()
  times.glmn     <- c()
  prec.glmn      <- c()
  tp.glmn <- c()
  tn.glmn <- c()
  fp.glmn <- c()
  fn.glmn <- c()

  for (thres in thres.glmn) {
##    cat(" ", thres)
    r.start <- proc.time()
    if (mine) {
      lasso.glmn <- lasso(x.train, y.train, lambda1=lambda, control=list(method="pathwise", threshold=thres, timer=TRUE))
    } else {
      lasso.glmn <- glm2fit(glmnet(x.train, y.train, lambda=lambda/sqrt(n),thresh=thres),n)
    }
    timer <- (proc.time() - r.start)[3]

    ##    cat(' ', length(lasso.glmn@lambda))
    error.glmn.mse <- cbind(error.glmn.mse, error(lasso.glmn, x.test, y.test))
    error.glmn.cla <- cbind(error.glmn.cla, error(lasso.glmn, x.test, y.test, type="class", thres=thres.class))
    error.glmn.sup <- cbind(error.glmn.sup,apply(lasso.glmn@coefficients, 1, function(x) sum(sign(x) != sign(true.beta))))

  tp.glmn <- cbind(tp.glmn, apply(lasso.glmn@coefficients, 1, function(x) length(intersect(which(x != 0), which(true.beta != 0)))))
  tn.glmn <- cbind(tn.glmn, apply(lasso.glmn@coefficients, 1, function(x) length(intersect(which(x == 0), which(true.beta == 0)))))
  fp.glmn <- cbind(fp.glmn, apply(lasso.glmn@coefficients, 1, function(x) length(intersect(which(x != 0), which(true.beta == 0)))))
  fn.glmn <- cbind(fn.glmn, apply(lasso.glmn@coefficients, 1, function(x) length(intersect(which(x == 0), which(true.beta != 0)))))
  
  objc <- objective(lasso.glmn,x.train,y.train)
    prec <- sqrt(sum( (objc - objc.ref)^2))/length(lambda)
    prec.glmn  <- c(prec.glmn , rep(prec, n.lambda))
    times.glmn <- c(times.glmn, rep(timer, n.lambda))
  }
  
##   cat("\nSPAMs for various thresholds...")
##   ## spams for various thresholds...
##   error.spam.mse <- c()
##   error.spam.cla <- c()
##   error.spam.sup <- c()
##   times.spam     <- c()
##   prec.spam      <- c()

##   for (thres in thres.spam) {
##     cat(" ", thres)
##     if (mine) {
##       lasso.spams <- elastic.net(x.train, y.train, lambda=lambda, gamma=gamma, control=list(method="fista", threshold=thres, timer=TRUE))
##     } else {
##       lasso.spams <- lasso.spam(x.train, y.train, lambda, type="fista", tol=thres.spam)
##     }
##     ##    cat(' ', length(lasso.spams@lambda))

##     error.spam.mse <- cbind(error.spam.mse, error(lasso.spams, x.test, y.test))
##     error.spam.cla <- cbind(error.spam.cla, error(lasso.spams, x.test, y.test, type="class", thres=thres.class))
##     error.spam.sup <- cbind(error.spam.sup, apply(lasso.spams@coefficients, 1, function(x) sum(sign(x) != sign(beta))))

##     objc <- objective(lasso.spams,x.train,y.train)
##     prec <- sqrt(sum( (objc - objc.ref)^2))/length(lambda)
##     prec.spam  <- c(prec.spam , rep(prec, n.lambda))
##     times.spam <- c(times.spam, rep(lasso.spams@monitoring$external.timer, n.lambda))
##   }

  error.mse  <- c(error.quad.mse, c(error.glmn.mse))
  error.cla  <- c(error.quad.cla, c(error.glmn.cla))
  error.sup  <- c(error.quad.sup, c(error.glmn.sup))
  tp <- c(tp.quad, c(tp.glmn))
  tn <- c(tn.quad, c(tn.glmn))
  fp <- c(fp.quad, c(fp.glmn))
  fn <- c(fn.quad, c(fn.glmn))

  times      <- c(times.quad, times.glmn)
  prec       <- c(prec.quad, prec.glmn)

  lambdas    <- rep(lambda, 1+length(thres.glmn))
  steps      <- rep(1:n.lambda, 1+length(thres.glmn))
  package    <- factor(rep(c("quad","glmn"), c(n.lambda, n.lambda*length(thres.glmn))), levels=c("quad","glmn"),ordered=TRUE)
  method     <- factor(rep(c("quad",paste("glmn",thres.glmn)), each=n.lambda), levels=c("quad",paste("glmn",thres.glmn)),ordered=TRUE)
  
  return(data.frame(steps=steps,lambda=lambdas,
                    error.mse=error.mse,error.cla=error.cla,error.sup=error.sup,tp=tp, tn=tn,fp=fp,fn=fn,
                    times=times,prec=prec, package=package,method=method))
} 

getData <- function(study = c("golub", "guedj", "lung", "julia", "prostate")) {

  studies <- c("golub", "guedj", "lung", "julia", "prostate")

  stopifnot(study %in% studies)

  if (study == "golub") {
    ## GOLUB DATA SET
    library("golubEsets")
    data(Golub_Test)
    data(Golub_Train)
    x <- rbind(t(exprs(Golub_Train)),
               t(exprs(Golub_Test)))
    y <- c(Golub_Train$ALL.AML,Golub_Test$ALL.AML)-1    
  }
  
  if (study == "guedj") {
    source("data/test.R")
    load("data/guedj11.RData")
    ## nombre de gènes choisis
    N <- 500
    out <- clearData(annot,guedj11)
    ann <- out$annot
    gue <- out$data
    y <- gue['205225_at',]
    ## y <- ann
    expr <- gue[!is.element(rownames(gue), c('205225_at', '211233_x_at', '211234_x_at','211235_s_at','211627_x_at','215551_at','215552_s_at','217163_at','217190_x_at')),]
    genes <- welch(ann,expr,crit.gene=N)
    x <- t(expr[genes[4,],])
  }

  if (study == "lung") {
    ## loading expression data
    load("datasets/expr_lung.RData")
    ## annotation file
    annot  <- read.delim("datasets/annotAffy.txt")
    ## loading phenotype data
    load("datasets/pheno_lung.RData")
    phenot <- rep(0, length(pheno))
    phenot[pheno == "tumor"] <- 1
    x  <- t(expr)
    y <- phenot
  }

  if (study == "julia") {
    ## DAS28 + Polyarthrite
    load("~/svn/select/data/Julia.RData")
    y <- dataExpressionSet$DAS28_wk0
    x <- t(exprs(dataExpressionSet))
    x <- x[, !is.na(colSums(x))]    
  }

  if (study == "prostate") {
    ## prostate data
    prostate <- read.table(file="http://www-stat.stanford.edu/~tibs/ElemStatLearn/datasets/prostate.data")
    x  <- as.matrix(prostate[, 1:8])
    y  <- c(prostate[, 9])
  }

  return(list(x=x,y=y))
}

get.mean.error <- function(out, alpha=0.05) {

  u <- qnorm(1-alpha/2)
  
  error.sup.mean  <- sapply(split(out, out$instance), function(x) tapply(x$error.sup, x$lambda, mean))
  error.sup.se    <- sapply(split(out, out$instance), function(x) tapply(x$error.sup, x$lambda, sd))
  error.sup.lower <- error.sup.mean - u*error.sup.se
  error.sup.upper <- error.sup.mean + u*error.sup.se

  error.mse.mean  <- sapply(split(out, out$instance), function(x) tapply(x$error.mse, x$lambda, mean))
  error.mse.se    <- sapply(split(out, out$instance), function(x) tapply(x$error.mse, x$lambda, sd))
  error.mse.lower <- error.mse.mean - u*error.mse.se
  error.mse.upper <- error.mse.mean + u*error.mse.se

  error.cla.mean  <- sapply(split(out, out$instance), function(x) tapply(x$error.cla, x$lambda, mean))
  error.cla.se    <- sapply(split(out, out$instance), function(x) tapply(x$error.cla, x$lambda, sd))
  error.cla.lower <- error.cla.mean - u*error.cla.se
  error.cla.upper <- error.cla.mean + u*error.cla.se

  lambda <- as.numeric(rownames(error.sup.mean))

  error.mean <- data.frame(lambda=lambda,
                           error.sup.mean=error.sup.mean, error.sup.upper=error.sup.upper, error.sup.lower=error.sup.lower,
                           error.mse.mean=error.mse.mean, error.mse.upper=error.mse.upper, error.mse.lower=error.mse.lower,
                           error.cla.mean=error.cla.mean, error.cla.upper=error.cla.upper, error.cla.lower=error.cla.lower)

  error.mean <- error.mean[order(lambda), ]

  return(error.mean)
}
