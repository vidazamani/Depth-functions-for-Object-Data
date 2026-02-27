//MOD2.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector MOD2_cpp(NumericMatrix D) {
  int n = D.nrow();
  NumericMatrix p(n, n);
  NumericVector ojadepth(n);
  NumericVector area(n);
  
  for (int k = 0; k < n; ++k) {
    // Reset p matrix for each k
    for (int i = 0; i < n - 1; ++i) {
      for (int j = i + 1; j < n; ++j) {
        double dik2 = D(i, k) * D(i, k);
        double djk2 = D(j, k) * D(j, k);
        double dij2 = D(i, j) * D(i, j);
        
        // Construct the 2x2 matrix S
        double a = dik2;
        double c = -0.5 * (dij2 - dik2 - djk2);
        double d = djk2;
        
        // Determinant of 2x2 matrix
        double det = a * d - c * c;
        
        // Filter small values
        if (det < 1e-6) det = 0.0;
        
        p(i, j) = det;
      }
    }
    
    // Sum sqrt of upper triangle of p
    double sum_sqrt = 0.0;
    for (int i = 0; i < n - 1; ++i) {
      for (int j = i + 1; j < n; ++j) {
        sum_sqrt += std::sqrt(p(i, j));
      }
    }
    
    area[k] = sum_sqrt;
    double norm_factor = 0.5 * n * (n - 1);
    ojadepth[k] = 1.0 / (1.0 + (area[k] / norm_factor));
  }
  
  return ojadepth;
}
