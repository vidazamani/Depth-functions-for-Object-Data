library(ICtest)
library(shapes)


# 5 x 5 covariance matrices
p <- 5

# A total of n = 100 random covariance matrices
n <- 100

many_S <- replicate(n, {
  
  U <- rorth(p)
  D <- diag(exp(rnorm(p)))
  S <- U%*%D%*%t(U)
  
  diagS <- diag(S)
  diag(1/sqrt(diagS))%*%S%*%diag(1/sqrt(diagS))
  
})


deepest_point <- diag(p)



# Distance matrix between the correlation matrices

library(CovTools)

D <- CovDist(many_S, method = "AIRM")






# Halfspace (Tukey) metric depth between training sample and test point
MHD_V2 <- function(D, d){
  
  n <- nrow(D)
  p = matrix(0,n,n)
  tu = c()
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      s = 0
      for (k in 1:n) {
        if (D[k,i]<= D[k,j]+1e-6){
          s=s+1
        }
      }
      
      p[i,j] = (1/n)*s
      
    }
  }
  

    Q = NULL
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        
        if (d[i]<=d[j]){
          Q = c(Q,p[i,j])
        }
      }
    }
    tu = min(Q)
 
  
  
  
  return(tu)
  
}





## Depth function for Test point

test_point_depth = function(v,D, many_S){
  
  vector_to_matrix = function(Cholvec){
    dim = sqrt(2*length(Cholvec)+(1/4)) - 0.5
    z = matrix(0,dim,dim)
    z[upper.tri(z, diag = TRUE)] <- Cholvec
    return(t(z)%*%z)
  }
  
  CM_test = vector_to_matrix(v) 
  
  d = c()
  for (i in 1:dim(many_S)[3]) {
    d[i] = distcov(many_S[,,i],CM_test)
  }
  
  MHD_V2(D,d)
}




Matrix_to_vector = function(CovMatrix){
  
  U <- chol(CovMatrix)
  U[upper.tri(U, diag = TRUE)]
  
}


testpoint = many_S[,,20]
v = Matrix_to_vector(testpoint)

test_point_depth(v,D, many_S)

