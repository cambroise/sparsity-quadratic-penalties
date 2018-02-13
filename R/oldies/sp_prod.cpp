// Rscript -e "RcppArmadillo:::SHLIB('sp_prod.cpp')"

#include <RcppArmadillo.h>

using namespace Rcpp ;

RcppExport SEXP sp_prod(SEXP X,
                        SEXP VAR,
                        SEXP N  ,
                        SEXP P) {

  Rcpp::List  SX       = Rcpp::List(X) ; // sparsely encoded SCALED design matrix
  arma::uword n        = Rcpp::as<int>       (N)         ; // sample size
  arma::uword p        = Rcpp::as<int>       (P)         ; // problem size
  arma::uword var      = Rcpp::as<int>       (VAR)       ; // column added
  arma::uvec Xi        ; // row indices of nonzeros
  arma::uvec Xj        ; // col indices of nonzeros
  arma::vec  Xx        ; // values of nonzeros
  Xi        = Rcpp::as<arma::uvec>(SX[0]);
  Xj        = Rcpp::as<arma::uvec>(SX[1]);
  Xx        = Rcpp::as<arma::vec> (SX[2]);
  
  arma::vec new_col = arma::zeros<arma::vec>(p) ;

  // non sparse encoding of the currently activated variable
  arma::vec var_col = arma::zeros<arma::vec>(n) ; 
  arma::uvec indXjvar = find(Xj == var);
  var_col.elem(Xi.elem(indXjvar)) = Xx.elem(indXjvar);

  arma::uvec indXjj ;
  for (int j=0; j<p; j++) {
    // non null entries of the current col j of X
    indXjj = find(Xj == j);
    new_col(j) = dot(Xx.elem(indXjj), var_col.elem(Xi.elem(indXjj)) );
  }
  
  return Rcpp::wrap(new_col) ;
}

RcppExport SEXP sp_prod2(SEXP X,
			 SEXP VAR,
			 SEXP N,
			 SEXP P) {

  Rcpp::List  SX       = Rcpp::List(X) ; // sparsely encoded SCALED design matrix
  arma::uword var      = Rcpp::as<int>       (VAR)       ; // column added
  arma::uword p        = Rcpp::as<int>       (P)         ; // problem size
  arma::uword n        = Rcpp::as<int>       (N)         ; // sample size
  arma::uvec Xi        ; // row indices of nonzeros
  arma::uvec Xp        ; // indices of first nonzero of each column
  arma::uvec Xnp       ; // # of nonzero in each column
  arma::vec  Xx        ; // values of nonzeros
  Xi        = Rcpp::as<arma::uvec>(SX[0]);
  Xp        = Rcpp::as<arma::uvec>(SX[1]);
  Xnp       = Rcpp::as<arma::uvec>(SX[2]);
  Xx        = Rcpp::as<arma::vec> (SX[3]);

  arma::vec new_col = arma::zeros<arma::vec>(p) ;

  // consider just the column of X which are non zero
  arma::uvec j_nz = find(Xnp > 0);
  arma:: vec col_Vx ;
  arma::uvec col_Xi ;
  arma:: vec col_Xx ;

  // if any nonzero in the X[, var] column, do the sparse product
  if (Xnp[var] > 0) {

    col_Vx = arma::zeros<arma::vec>(n) ;
    col_Vx.elem(Xi.subvec(Xp[var],Xp[var+1]-1)) = Xx.subvec(Xp[var],Xp[var+1]-1);
    
    // loop along each column of X
    for (int j=0; j<j_nz.n_elem; j++) {
      col_Xx = Xx.subvec(Xp[j_nz[j]],Xp[j_nz[j]+1]-1) ;
      col_Xi = Xi.subvec(Xp[j_nz[j]],Xp[j_nz[j]+1]-1) ;
      new_col[j_nz[j]] = dot(col_Xx, col_Vx.elem(col_Xi));
    }
  }
  
  return Rcpp::wrap(new_col) ;
}
