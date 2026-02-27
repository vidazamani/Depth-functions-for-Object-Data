//MLD.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector MLD_cpp(NumericMatrix D) {
  int n = D.nrow();
  
  NumericVector lens_depth(n);
  
  for (int p = 0; p < n; ++p) {
    int s = 0;

    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        double maximum = std::max(D(i, p), D(j, p));
        if (D(i, j) > maximum + 1e-6) {
          ++s;
        }
      }
    }

    lens_depth[p] = s / (n * (n - 1) / 2.0);
  }
  
  return(lens_depth);
}
