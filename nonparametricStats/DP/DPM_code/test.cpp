// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <Rmath.h>
using namespace Rcpp; using namespace arma;


// [[Rcpp::export]]
arma::vec test(arma::vec x, int i) {
  int N = x.size();
  arma::vec temp = x.head(i-1);
  arma::vec temp2 = x.tail(N-i);
  int i1 = temp.size();
  int i2 = temp2.size();
  int i3 = i1 + i2;
  arma::vec temp3(i3);
  for (int j = 0; j < i1; j++) {
    temp3(j) = temp(j);
  }

  for (int k = 0; k < i2; k++) {
    temp3(k + i1) = temp2(k);
  }
  return temp3;
}



arma::uvec test2(arma::uvec x, int i) {
  int N = x.size();
  arma::uvec temp = x.head(i-1);
  arma::uvec temp2 = x.tail(N-i);
  int i1 = temp.size();
  int i2 = temp2.size();
  int i3 = i1 + i2;
  arma::uvec temp3(i3);
  for (int j = 0; j < i1; j++) {
    temp3(j) = temp(j);
  }

  for (int k = 0; k < i2; k++) {
    temp3(k + i1) = temp2(k);
  }
  return temp3;
}
