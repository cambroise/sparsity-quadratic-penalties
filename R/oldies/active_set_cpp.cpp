// Rscript -e "RcppArmadillo:::SHLIB('active_set_cpp.cpp')"

#include <RcppArmadillo.h>
#include "quadra_cpp.h"
#include "fista_cpp.h"
#include "pathwise_cpp.h"

RcppExport SEXP active_set_cpp(SEXP X0, SEXP XTX, SEXP X, SEXP D, SEXP XTY, SEXP GRD, SEXP LAMBDA, SEXP EPS, SEXP FUN) ;

using namespace Rcpp ;

RcppExport SEXP active_set_cpp(SEXP X0, SEXP XTX, SEXP X, SEXP D, SEXP XTY, SEXP GRD, SEXP LAMBDA, SEXP EPS, SEXP FUN) {

  // Reading input variables
  arma::vec xk  = Rcpp::as<arma::vec>(X0)      ; // vector of parameters
  arma::mat xtx                                ; // covariance matrix
  arma::mat xtx_A                              ; // covariance matrix of the activated variable
  arma::mat d   = Rcpp::as<arma::mat>(D)       ; // ridge matrix
  arma::mat x   = Rcpp::as<arma::mat>(X)       ; // design matrix
  if (!xk.is_empty()) {
    arma::mat xtx = Rcpp::as<arma::mat>(XTX)   ; 
  }
  arma::vec xty   = Rcpp::as<arma::vec>(XTY)   ; // reponses to predictors vector
  arma::vec grd   = Rcpp::as<arma::vec>(GRD)   ; // smooth part of the gradient
  arma::uvec A    = find(xk)                   ; // set of currently activated variables
  double lambda   = Rcpp::as<double>(LAMBDA)   ; // penalty level
  double eps      = Rcpp::as<double>(EPS)      ; // precision required
  int    fun      = Rcpp::as<int>(FUN)         ; // the chosen solver (1=pathwise,2=proximal,0=quadra)

  std::cout << "data read \n";
      
  // Initializing variables 
  int iter     = 0      ; // current iterate
  int max_iter = 5000   ; // max # of iterates
  int var_in            ; // currently added variable
  int nbr_in = 0        ; // # of currently added variables
  arma::uvec is_in      ; // a variable to check if a variable is already in the active set
  Rcpp::List out_optim  ; // the list of output of the optimization function
  arma::uvec iter_optim ; // a vector which stores the # of iterates along the optimization 
  arma::vec x1          ; // first  solution of the optimization routine
  arma::vec x2          ; // second solution of the optimization routine (only for quadra)
  arma::vec grd1        ; // gradient of the first solution 
  arma::vec grd2        ; // gradient of the second solution 
  arma::vec gaps_swap   ; // dual gap for unactive variable
  arma::uvec swap       ; // stores the variables whose sign has swaped during optimization
  arma::uvec null       ; // stores the variables which go to zero during optimization

  std::cout << "variable initialized \n";
  
  arma::vec gaps = abs(grd) - lambda                        ; // dual gap for unactive variable
  arma::uvec gaps_null    = find(gaps < 0)                  ;
  if (gaps_null.n_elem > 0) {gaps.elem(gaps_null).fill(0);} ;
  gaps.elem(A) = abs(grd.elem(A) + lambda * xk/abs(xk))     ; // dual gap for active variable
  double gap = arma::as_scalar(max(gaps))                   ; // dual gap of the problem
  double gap_old                                            ;
  
  std::cout << "initial gap values computed \n";
  
  // monitoring algorithm convergence (dual gap condition, max # of iterate reached, gap stability)
  Rcpp::LogicalVector conv = LogicalVector::create(gap<eps, 0, 0); ; 
   
  while (!is_true(any(conv))) {
    // _____________________________________________________________
    //
    // (1) VARIABLE ACTIVATION IF APPLICABLE
    // _____________________________________________________________
    var_in = as_scalar(find(gaps == max(gaps))) ;
    
    // Check if the variable is already in the active set
    is_in = find(A == var_in);
    if (is_in.n_elem == 0) {
      // If this is a newly added variable, then
      A.reshape(nbr_in+1,1)  ; // update the active set
      A[nbr_in] = var_in     ;
      xk.reshape(nbr_in+1,1) ; // update the vector of active parameters
      xk[nbr_in]    = 0.0    ;

      // Update the xtx and xtx_A matrices
      if (nbr_in > 0) {
	xtx_A = join_cols(xtx_A, xtx.row(var_in)) ;
      }
      xtx   = join_rows(xtx, strans(x) * x.col(var_in) + d.col(var_in)) ;
      xtx_A = join_rows(xtx_A, strans(xtx.row(var_in))) ;
      
      std::cout << "newly added variable " << var_in << "\n";
      nbr_in++;
    } else {
      std::cout << var_in << " is already in the active set - keep on optimizing" << "\n";
    }
    
    // _____________________________________________________________
    //
    // (2) OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
    // _____________________________________________________________
    //

    // std::cout << "Optimizing on\n";
    // std::cout << "A =" << strans(A) << "\n";
    //    std::cout << "xk =" << xk << "\n";
    //    std::cout << "xtx_A =" << xtx_A << "\n";
    switch (fun) {
    case 1 :
      out_optim = pathwise_cpp(xk, xtx_A, xty.elem(A), lambda, eps); 
      break;
    case 2 :
      out_optim = fista_cpp(xk, xtx_A, xty.elem(A), lambda, eps); 
      break;
    default: 
      out_optim = quadra_cpp(xk, xtx_A, xty.elem(A), lambda, eps);
    }
    // save the number of iterates performed along the optimization process
    iter_optim.reshape(iter+1,1) ;
    iter_optim[iter] = as_scalar(as<arma::uvec>(out_optim[0])) + 1 ;

    // _____________________________________________________________
    //
    // (3) VARIABLE DELETION IF APPLICABLE
    // _____________________________________________________________
    //
    if (out_optim.size() > 2) {
      x1 = as<arma::vec>(out_optim[1]);
      x2 = as<arma::vec>(out_optim[2]);
      swap = as<arma::uvec>(out_optim[3]);
      grd1 = -xty + xtx * x1;
      grd2 = -xty + xtx * x2;
      gaps_swap = abs(grd2.elem(A.elem(swap))  + lambda * x2.elem(swap)/abs(x2.elem(swap)));
      if (min(gaps_swap) < eps) {
	xk  = x2;
	grd = grd2;
      } else {
	xk  = x1;
	grd = grd1;
	std::cout << "removing variable" << swap << "\n";
	for (int j=0; j<swap.n_elem; j++) {
	  A.shed_row(swap[j]);
	  xk.shed_row(swap[j]) ;
	  xtx.shed_col(swap[j]) ;
	  xtx_A.shed_col(swap[j]) ;
	  xtx_A.shed_row(swap[j]) ;
	  nbr_in--;
	}
      }
    } else {
      x1  = as<arma::vec>(out_optim[1]);
      xk  = x1;
      grd = -xty + xtx * xk;
      null = find(abs(xk) + (abs(grd.elem(A)) - lambda) < 2*eps) ;
      if (!null.is_empty()) {
	std::cout << "removing variable" << null << "\n";
	for (int j=0; j<null.n_elem; j++) {
	  A.shed_row(null[j]);
	  xk.shed_row(null[j]) ;
	  xtx.shed_col(null[j]) ;
	  xtx_A.shed_col(null[j]) ;
	  xtx_A.shed_row(null[j]) ;
	  nbr_in--;
	}
      }
    }
    // std::cout << "done with xk =" << xk << "\n";

    // _____________________________________________________________
    //
    // (4) OPTIMALITY TESTING
    // _____________________________________________________________
    
    gaps         = abs(grd) - lambda                          ; // dual gap for unactive variables
    gaps_null    = find(gaps < 0)                             ;
    if (gaps_null.n_elem > 0) {gaps.elem(gaps_null).fill(0);} ;
    gaps.elem(A) = abs(grd.elem(A) + lambda * xk/abs(xk))     ; // dual gap for active variables
    gap_old      = gap                                        ;
    gap          = max(gaps)                                  ; // dual gap of the problem

    // Moving to the next iterate
    iter++;
    conv[0] = gap < eps; 
    conv[1] = iter > max_iter; 
    //conv[2] = abs(gap_old - gap) < eps; 
    conv[2] = 0;

    R_CheckUserInterrupt();    
  }

  return Rcpp::List::create(Rcpp::Named("beta.A")  = xk  ,
 			    Rcpp::Named("xtx.A")   = xtx ,
 			    Rcpp::Named("nabla.f") = grd ,
 			    Rcpp::Named("active")  = A   ,
 			    Rcpp::Named("gap")     = gap ,
 			    Rcpp::Named("conv")    = conv,
			    Rcpp::Named("iter")    = iter_optim);
    
}
