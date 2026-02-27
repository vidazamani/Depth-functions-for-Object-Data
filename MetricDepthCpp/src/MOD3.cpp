//MOD3.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double det3times3_cpp(NumericMatrix M) {
  return(M(0, 0) * (M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1)) -
         M(0, 1) * (M(1, 0) * M(2, 2) - M(1, 2) * M(2, 0)) +
         M(0, 2) * (M(1, 0) * M(2, 1) - M(1, 1) * M(2, 0)));
}

// [[Rcpp::export]]
NumericVector MOD3_cpp(NumericMatrix D) {
  int n = D.ncol();
  NumericVector ojadepth(n);
  NumericVector area(n);
  
  for (int w = 0; w < n; ++w) {
    double sum_sqrt = 0.0;
    
    for (int i = 0; i < (n - 2); ++i) {
      for (int j = i + 1; j < (n - 1); ++j) {
        for (int k = j + 1; k < n; ++k) {
          double diw2 = D(i, w) * D(i, w);
          double djw2 = D(j, w) * D(j, w);
          double dkw2 = D(k ,w) * D(k, w);
          double dij2 = D(i, j) * D(i, j);
          double dik2 = D(i, k) * D(i, k);
          double djk2 = D(j, k) * D(j, k);
          
          NumericVector M = {
            diw2,
             -0.5 * (dij2 - djw2 - diw2),
             -0.5 * (dik2 - dkw2 - diw2),
             -0.5 * (dij2 - diw2 - djw2),
              djw2,
              -0.5 * (djk2 - dkw2 - djw2),
              -0.5 * (dik2 - diw2 - dkw2),
               -0.5 * (djk2 - djw2 - dkw2),
               dkw2
          };
          
          M.attr("dim") = Dimension(3, 3);
          
          double detM = M(0, 0) * (M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1)) -
            M(0, 1) * (M(1, 0) * M(2, 2) - M(1, 2) * M(2, 0)) +
            M(0, 2) * (M(1, 0) * M(2, 1) - M(1, 1) * M(2, 0));
          sum_sqrt += std::sqrt(detM + 4.0 * diw2 * djw2 * dkw2);
        }
      }
    }
    
    area[w] = sum_sqrt;
    double norm_factor = (n - 2) * (n - 1) * n / 6;
    ojadepth[w] = 1.0 / (1.0 + (area[w] / norm_factor));
  }
  
  return ojadepth;
}


// [[Rcpp::export]]
double MOD3_test_cpp(NumericMatrix D, NumericVector d) {
  int n = D.ncol();
  double ojadepth = 0.0;
  double area = 0.0;
  
  d = d * d;
  

  double sum_sqrt = 0.0;
    
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        double dij2 = D(i, j) * D(i, j);
        double dik2 = D(i, k) * D(i, k);
        double djk2 = D(j, k) * D(j, k);
        
        NumericVector M = {
          d(i),
          -0.5 * (dij2 - d(j) - d(i)),
          -0.5 * (dik2 - d(k) - d(i)),
          -0.5 * (dij2 - d(i) - d(j)),
          d(j),
          -0.5 * (djk2 - d(k) - d(j)),
          -0.5 * (dik2 - d(i) - d(k)),
          -0.5 * (djk2 - d(j) - d(k)),
          d(k)
        };
        
        M.attr("dim") = Dimension(3, 3);
        
        double detM = M(0, 0) * (M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1)) -
          M(0, 1) * (M(1, 0) * M(2, 2) - M(1, 2) * M(2, 0)) +
          M(0, 2) * (M(1, 0) * M(2, 1) - M(1, 1) * M(2, 0));
        sum_sqrt += std::sqrt(detM + 4.0 * d(i) * d(j) * d(k));
      }
    }
  }
  
  area = sum_sqrt;

  double norm_factor = 0.5 * (n - 1) * n * n;
  ojadepth = 1.0 / (1.0 + (area / norm_factor));

  return(ojadepth);
}

