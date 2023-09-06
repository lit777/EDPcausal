// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dot_prod(NumericMatrix xmat, 
                       NumericMatrix beta_mat, 
                       int n, int p) {
  
  NumericVector res(n);
  NumericVector xvec(p);
  NumericVector beta_vec(p);
  
  for(int i=0; i < n; ++i){
    xvec = xmat(i, _ );
    beta_vec = beta_mat(_ , i);
  
    res(i) = sum(xvec*beta_vec);
  }
  
  return res;
}