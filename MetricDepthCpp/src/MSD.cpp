//MSD.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector MSD_cpp(NumericMatrix D) {
  int n = D.ncol();
  NumericVector res(n);
  
  for (int k = 0; k < n; ++k) {
    double res_now = 0.0;
    NumericVector d_row = D.row(k);
    
    // Create index_set = {0, ..., n-1} \ {k} \ {i | d_row[i] < 1e-6}
    NumericVector index_set;
    for (int i = 0; i < n; ++i) {
      if (i != k && d_row[i] >= 1e-6) {
        index_set.push_back(i);
      }
    }
    
    for (int i : index_set) {
      for (int j : index_set) {
        double temp = d_row[i] / d_row[j];
        res_now += temp + 1.0 / temp - std::pow(D(i ,j), 2) / (d_row[i] * d_row[j]);
      }
    }
    
    res[k] = res_now;
  }
  
  // Normalize and transform result
  for (int k = 0; k < n; ++k) {
    res[k] = 1.0 - 0.5 * res[k] / (n * n);
  }
  
  return res;
}
