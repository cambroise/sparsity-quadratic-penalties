rm(list=ls())


library(RcppArmadillo)
library(inline)

src <-'
arma::mat R = Rcpp::as<arma::mat>(Rs) ;
arma::mat X = Rcpp::as<arma::mat>(Xs) ;

int p = X.n_cols;
arma::colvec rp = arma::zeros<arma::colvec>(p,1);
rp.rows(0,p-2) = solve (trimatl(strans(R)), strans(X.cols(0,p-2)) * X.col(p-1) );
rp[p-1] = sqrt(dot(X.col(p-1),X.col(p-1))  - dot(rp,rp));
R = join_rows( join_cols(R, arma::zeros<arma::mat>(1,p-1)) , rp);

return(wrap(R));
'

cholupdate.arma <- cxxfunction(signature(Rs="matrix", Xs="matrix"),
                               src,
                               plugin="RcppArmadillo")

cholupdate <- function (R, X) {
  ## Update a Cholesky decomposition
  if (is.null(R)) {
    R <- matrix(sqrt(sum(X^2)), 1, 1)
  } else {
    p <- ncol(X)
    r.p <- backsolve(R, crossprod(X[,1:(p-1)], X[, p]), transpose=TRUE)
    rpp <- sqrt(sum(X[,p]^2) - sum(r.p^2))
    R <- cbind(rbind(R, 0), c(r.p, rpp))
  }
  return(R)
}

src <-'
arma::mat R = Rcpp::as<arma::mat>(Rs) ;
int j = Rcpp::as<int>(js) -1;
double eps = 1e-10;

arma::vec x = arma::zeros<arma::vec>(2,1);
arma::mat G = arma::zeros<arma::mat>(2,2);

R.shed_col(j);
int p = R.n_cols;
double r;
for (int k=j; k<p; k++) {
   x = R.submat(k,k,k+1,k);

   if (x[1] != 0) {
     r = norm(x,2);
     G <<  x(0) << x(1) << arma::endr
       << -x(1) << x(0) << arma::endr;
     G = G / r;
     x(0) = r; x(1) = 0;
   } else {
     G = arma::eye(2,2);
   }
   R.submat(k,k,k+1,k) = x;
   if (k < p-1) {
      R.submat(k,k+1,k+1,p-1) = G * R.submat(k,k+1,k+1,p-1);
   }
  }
  R.shed_row(p);

return(wrap(R));
'
choldowndate.arma <- cxxfunction(signature(Rs="matrix", js="numeric"),
                                 src,
                                 plugin="RcppArmadillo")

planerot <- function(x) {
  if (x[2] != 0) {
    r <- sqrt(sum(x^2))
    G <- matrix(c(x,-x[2],x[1]),2,2, byrow=TRUE)/r
    x <- c(r,0)
  } else {
    G <- matrix(c(1,0,0,1),2,2)
  }
  return(list(G=G,x=x))
}

choldowndate <- function (R, j) {
  ## Downdate a Cholesky decomposition
  R <- R[,-j]
  p <- ncol(R)
  for (k in j:p) {
    n <- k:(k+1)
    out <- planerot(R[n,k])
    R[n,k] <- out$x
    if (k < p)
      R[n,(k+1):p] <- out$G %*% R[n,(k+1):p]
  }
  R <- R[-(p+1),]
  return(R)
}

cholupdateGram <- function (R, XtX) {
  ## Update a Cholesky decomposition
  if (is.null(R)) {
    R <- matrix(sqrt(XtX), 1, 1)
  } else {
    p <- ncol(XtX)
    r.p <- backsolve(R, XtX[-p,p], transpose=TRUE)
    rpp <- sqrt(XtX[p,p] - sum(r.p^2))
    R <- cbind(rbind(R, 0), c(r.p, rpp))
  }
  return(R)
}

src <-'
arma::mat R   = Rcpp::as<arma::mat>(Rs) ;
arma::mat XtX = Rcpp::as<arma::mat>(XtXs) ;

int p = XtX.n_cols;

arma::colvec rp  = arma::zeros<arma::colvec>(p,1);
rp.rows(0,p-2) = solve (trimatl(strans(R)), XtX.submat(0,p-1,p-2,p-1));
rp(p-1) = sqrt(XtX(p-1,p-1) - dot(rp,rp));
R = join_rows( join_cols(R, arma::zeros<arma::mat>(1,p-1)) , rp);

return(wrap(R));
'

cholupdate.arma2 <- cxxfunction(signature(Rs="matrix", XtXs="matrix"),
                               src,
                               plugin="RcppArmadillo")

p <- 200
n <- 500
x <- matrix(rnorm(p*n),n,p)
x <- scale(x)
xtx <- crossprod(x,x)

time.chol <- system.time(for (i in 1:p) {
  chol.R <- chol(t(x[,1:i]) %*% x[,1:i])
})

chol.R.X <- NULL
time.R.X <- system.time(for (i in 1:p) {
  chol.R.X <- cholupdate(chol.R.X,x[,1:i])
})

chol.R.XtX <- NULL
time.R.XtX <- system.time(for (i in 1:2) {
  chol.R.XtX <- cholupdateGram(chol.R.XtX,xtx[1:i,1:i])
})

chol.arma.X <- matrix(sqrt(sum(x[,1]^2)),1,1)
time.arma.X <- system.time(for (i in 2:p) {
  chol.arma.X <- cholupdate.arma(chol.arma.X,x[,1:i])
})

chol.arma.XtX <- matrix(sqrt(xtx[1,1]),1,1)
time.arma.XtX <- system.time(for (i in 2:p) {
  chol.arma.XtX <- cholupdate.arma2(chol.arma.XtX,xtx[1:i,1:i])
})


chol.R.dn <- choldowndate(chol.R, 57)
chol.arma.XtX.dn <- choldowndate.arma(chol.arma.XtX, 57)

D.dn <- sum((chol.arma.XtX.dn - chol.R.dn)^2)

D.up <- sum((chol.arma.X - chol.R)^2)
D.up2 <- sum((chol.arma.XtX - chol.R)^2)
