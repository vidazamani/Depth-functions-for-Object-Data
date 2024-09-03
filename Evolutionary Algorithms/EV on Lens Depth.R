library(ICtest)
library(shapes)

set.seed(1)

# 5 x 5 corr matrices
p <- 5

# A total of n = 100 random corr matrices
n <- 10

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





# Lens depth for Object data between train sample and test point

EMD_v3 <- function(D, d){
  
  n <- nrow(D)
  lens_depth = c()
  s = 0
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
          
          if (D[i,j]> max(d[i],d[j])+1e-6) {
            s = s+1
          }
      }
      
    }
    
    
    lens_depth  = (1/choose(n,2))*s
    
  
  
  return(lens_depth)
  
}



EMD <- function(D){
  
  n <- nrow(D)
  lens_depth = c()
  
  for (k in 1:n) {
    s = 0
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        
          if (D[i,j]> max(D[i,k],D[j,k])+1e-6) {
            s = s+1
           
          }
          
      
        
      }
      
    }
    
    
    lens_depth[k]  = (1/choose(n,2))*s
    
  }
  
  return(lens_depth)
  
}



## Depth function for test point

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
  
  EMD_v3(D,d)
}




# D[,6]
# d
# 
# EMD(D)
# 
# EMD_v3(D,d)
# 




Matrix_to_vector = function(CovMatrix){
  
  U <- chol(CovMatrix)
  U[upper.tri(U, diag = TRUE)]
  
}

testpoint = many_S[,,8]
v = Matrix_to_vector(testpoint)

round(test_point_depth(v,D, many_S),3)


round(EMD(D),3)

