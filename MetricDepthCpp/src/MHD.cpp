//MHD.cpp
#include <Rcpp.h>
using namespace Rcpp;





// [[Rcpp::export]]
NumericVector MHD_cpp(NumericMatrix D) {
  int n = D.nrow();
  NumericMatrix p(n, n);
  NumericVector tu(n);
  
  // Compute matrix p
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double s = 0.0;
      for (int k = 0; k < n; ++k) {
        if (D(k, i) <= D(k, j) + 1e-6) {
          ++s;
        }
      }
      p(i, j) = s / n;
    }
  }
  
  // Compute tu vector
  for (int y = 0; y < n; ++y) {
    double Q = 1.0;
    for (int i = 0; i < n - 1; ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (D(y, i) <= D(y, j)) {
          Q = std::min(Q, p(i, j));
        }
      }
    }
    tu[y] = Q;
  }
  
  return tu;
}




// [[Rcpp::export]]
double MHD_test_cpp(NumericMatrix D, NumericVector d) {
  int n = D.ncol();
  NumericMatrix p(n, n);
  double Q = 1.0;
  
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double s = 0.0;
      for (int k = 0; k < n; ++k) {
        if (D(k,i) <= D(k,j) + 1e-6) {
          ++s;
        }
      }
      p(i,j) = s/n;
    }
  }
  
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (d[i] <= d[j]) {
        if (p(i,j) <= Q) {
          Q = p(i,j);
        }
      }
    }
  }
  
  return Q;
}
