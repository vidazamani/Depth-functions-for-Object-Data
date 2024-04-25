
## Loading necessary libraries 

library(ICtest)
library(abind)
library(CovTools)
library(ggplot2)
library(dplyr)
library(shapes)



# 5 x 5 covariance matrices
p <- 5

# A total of n = 100 random covariance matrices
n <- 50

# Proportion of outliers
eps <- 0.05

# Non-outlying sample size
n0 <- floor((1 - eps)*n)

# Outlying sample size
n1 <- n - n0

# How outlying the outliers are
mu <- 3


many_S_non_outliers <- replicate(n0, {
  
  U <- rorth(p)
  D <- diag(exp(rnorm(p)))
  S <- U%*%D%*%t(U)
  
  # diagS <- diag(S)
  # diag(1/sqrt(diagS))%*%S%*%diag(1/sqrt(diagS))
  # 
  S
})


many_S_outliers <- replicate(n1, {
  
  U <- rorth(p)
  D <- diag(exp(rnorm(p, mu)))
  S <- U%*%D%*%t(U)
  # 
  # diagS <- diag(S)
  # diag(1/sqrt(diagS))%*%S%*%diag(1/sqrt(diagS))
  # 
  S
})


many_S <- abind(many_S_non_outliers, many_S_outliers)
many_S[,,1]

# What is the deepest point of the non-outlying part of the sample on the population level?
# Because the "typical" value of N(0, 1) is 0, this makes the most typical value of D to be the identity matrix
# Consequently, the typical values of S and corr are also I (?)



deepest_point <- diag(p)




# 2. Compute the distance matrix
# Distance matrix between the correlation matrices


D <- CovDist(many_S, method = "AIRM")





# Metric three-dimensional Oja depth

Oja_v2 <- function(D){
  
  n <- nrow(D)
  B = array(0,c(n,n,n))
  ojadepth = c()
  area = c(0)
  
  
  for (w in 1:n) {
    
    for (i in 1:(n-1)) {
      for (j in i:n) {
        for (k in 1:n) {
          
          S = matrix(c((D[i,w])^2,
                       -0.5*((D[j,i])^2 - (D[j,w]^2) - (D[i,w])^2),
                       -0.5*((D[k,i])^2 - (D[k,w]^2) - (D[i,w])^2),
                       -0.5*((D[i,j])^2 - (D[i,w]^2) - (D[j,w])^2),
                       (D[j,w])^2,
                       -0.5*((D[k,j])^2 - (D[k,w]^2) - (D[j,w])^2),
                       -0.5*((D[i,k])^2 - (D[i,w]^2) - (D[k,w])^2),
                       -0.5*((D[j,k])^2 - (D[j,w]^2) - (D[k,w])^2),
                       (D[k,w])^2 )
                     ,3,3)
          
          B[i,j,k] = det(S)
          
        }
      }
    }
    B[B<0] = 0
    area[w] = sum(sqrt(B))
    ojadepth[w] = 1 / (1 + (1/(0.5*(n-1)*(n^2))*(area[w])))
  }
  
  return(ojadepth)
  
}



Oja_v2(D)

## Metric three-dimensional Oja depth between train sample and test point
Oja_v3 <- function(D, d){
  
  n <- nrow(D)
  B = array(0,c(n,n,n))
  ojadepth = c()
  area = c(0)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      for (k in 1:n) {
        
        S = matrix(c((d[i])^2,
                     -0.5*((D[j,i])^2 - (d[j]^2) - (d[i])^2),
                     -0.5*((D[k,i])^2 - (d[k]^2) - (d[i])^2),
                     -0.5*((D[i,j])^2 - (d[i]^2) - (d[j])^2),
                     (d[j])^2,
                     -0.5*((D[k,j])^2 - (d[k]^2) - (d[j])^2),
                     -0.5*((D[i,k])^2 - (d[i]^2) - (d[k])^2),
                     -0.5*((D[j,k])^2 - (d[j]^2) - (d[k])^2),
                     (d[k])^2 )
                   ,3,3)
        
        B[i,j,k] = det(S)
        
      }
    }
  }
  B[B<0] = 0
  area = sum(sqrt(B))
  ojadepth = 1 / (1 + (1/(0.5*(n-1)*(n^2))*(area)))
  
  
  return(ojadepth)
  
}




## Fitness function

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
  
  Oja_v3(D,d)
}




Matrix_to_vector = function(CovMatrix){
  
  U <- chol(CovMatrix)
  U[upper.tri(U, diag = TRUE)]
  
}


testpoint = many_S[,,50]
v = Matrix_to_vector(testpoint)

round(test_point_depth(v,D, many_S),3)








