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


## Oja depth between train sample and test point
Oja_v3 <- function(D, d){
  
  n <- nrow(D)
  p = matrix(0,n,n)
  ojadepth = c()
  area = c(0)
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        S = matrix(c((d[i])^2,
                     -0.5*((D[j,i])^2 - (d[j]^2) - (d[i])^2),
                     -0.5*((D[i,j])^2 - (d[i]^2) - (d[j])^2),
                     (d[j])^2)
                   ,2,2)
        
        p[i,j] = det(S)
        
      }
    }
    p[p<0] = 0 ### max(0,B) in paper
    area = sum(sqrt(p))
    ojadepth = 1 / (1 + (1/(0.5*(n^2-n))*(area)))
 
  
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


testpoint = many_S[,,10]
v = Matrix_to_vector(testpoint)

round(test_point_depth(v,D, many_S),3)

# round(Oja_v2(D),3)

set.seed(12)
ord = order(Oja_v2(D))


library(GA)

suggestedSol <- rbind(Matrix_to_vector(many_S[,,ord[n]]),
                      Matrix_to_vector(many_S[,,ord[n-1]]),
                      Matrix_to_vector(many_S[,,ord[n-2]]),
                      Matrix_to_vector(many_S[,,ord[n-3]]),
                      Matrix_to_vector(many_S[,,ord[n-4]]))



lower_1 = Matrix_to_vector(many_S[,,ord[n]]) - 0.5
upper_1 = Matrix_to_vector(many_S[,,ord[n]]) + 0.5

GA1 <- ga(type = "real-valued", 
          fitness =  function(x) fitness(x,D, many_S),
          lower = lower_1, 
          upper = upper_1, 
          suggestions = suggestedSol,
          popSize = 50, maxiter = 30)










GA1 <- ga(type = "real-valued", 
          fitness =  function(x) fitness(x,D, many_S),
          lower = c(rep(min(Matrix_to_vector(many_S[,,ord[n]])),0.5*(p*(p+1)))), 
          upper = c(rep(1.5,0.5*(p*(p+1)))), 
          suggestions = suggestedSol,
          popSize = 50, maxiter = 30)
summary(GA1)
plot(GA1)




GA1@solution

GA_So = vector_to_matrix(GA1@solution)


diagS <- diag(GA_So)
GA_Corr = diag(1/sqrt(diagS))%*%GA_So%*%diag(1/sqrt(diagS))


Diff <- array(0, dim = c(p, p, 2))


Diff[, , 1] <- GA_Corr
Diff[, , 2] <- diag(p)

out_of_sample_error = CovDist(Diff, method = "AIRM")[1,2]


Diff <- array(0, dim = c(p, p, 2))


Diff[, , 1] <- many_S[,,ord[n]]
Diff[, , 2] <- diag(p)

in_sample_error = CovDist(Diff, method = "AIRM")[1,2] ## this is exactly the "ojav2diff1" that we had


