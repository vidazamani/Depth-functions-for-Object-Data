library(ICtest)
library(shapes)
library(factoextra)



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




# Metric Spatial depth between train sample and test point

MSD_V3 <- function(D, d){
  
  n <- nrow(D)

    
    res_now <- 0
    
    d_row <- d
    
    index_set <- setdiff(1:n, which(d_row < 1e-6))
    
    for(i in index_set){
      
      for(j in index_set){
        
        temp <- d_row[i]/d_row[j]
        
        res_now <- res_now + temp + 1/temp - D[i, j]^2/(d_row[i]*d_row[j])
        
      }
      
    }
    

  
    res_now <- res_now/n^2
  
  1 - 0.5*res_now
  
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
  
  MSD_V3(D,d)
}




Matrix_to_vector = function(CovMatrix){
  
  U <- chol(CovMatrix)
  U[upper.tri(U, diag = TRUE)]
  
}


testpoint = many_S[,,10]
v = Matrix_to_vector(testpoint)

test_point_depth(v,D, many_S)


MSD(D)



