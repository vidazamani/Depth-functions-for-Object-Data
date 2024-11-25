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



Matrix_to_vector = function(CovMatrix){
  
  U <- chol(CovMatrix)
  U[upper.tri(U, diag = TRUE)]
  
}


testpoint = many_S[,,2]
v = Matrix_to_vector(testpoint)


#### PCA
encoded_train = matrix(0,dim(many_S)[3],length(v))

for (i in 1:dim(many_S)[3]) {
  encoded_train[i,] = Matrix_to_vector(many_S[,,i])
}
res.pca <- prcomp(encoded_train, scale = FALSE)
fviz_eig(res.pca, ncp = 15)


loading = res.pca$rotation

pca_testpoint = as.vector(t(loading)%*%(v-colMeans(encoded_train)))


#### PCA Fitness function


test_point_depth = function(v,D, many_S){
  
  vector_to_matrix = function(Cholvec){
    dim = sqrt(2*length(Cholvec)+(1/4)) - 0.5
    z = matrix(0,dim,dim)
    z[upper.tri(z, diag = TRUE)] <- Cholvec
    return(t(z)%*%z)
  }
  
  v_star = as.vector(loading[,1:length(v)]%*%v + colMeans(encoded_train))
                              
  
  CM_test = vector_to_matrix(v_star) 
  
  d = c()
  for (i in 1:dim(many_S)[3]) {
    d[i] = distcov(many_S[,,i],CM_test)
  }
  
  Oja_v3(D,d)
}

test_point_depth(pca_testpoint[1:4],D, many_S)


# Oja depth 2nd Version
Oja_v2 <- function(D){
  
  n <- nrow(D)
  p = matrix(0,n,n)
  ojadepth = c()
  area = c(0)
  
  
  for (k in 1:n) {
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        S = matrix(c((D[i,k])^2,
                     -0.5*((D[j,i])^2 - (D[j,k]^2) - (D[i,k])^2),
                     -0.5*((D[i,j])^2 - (D[i,k]^2) - (D[j,k])^2),
                     (D[j,k])^2)
                   ,2,2)
        
        p[i,j] = det(S)
        
      }
    }
    #area[k] = sum(upper.tri(p)*p)
    area[k] = sum(sqrt(p))
    ojadepth[k] = 1 / (1 + (1/(0.5*(n^2-n))*(area[k])))
  }
  
  return(ojadepth)
  
}

library(GA)

set.seed(12)
ord = order(Oja_v2(D))

suggestedSol <- rbind(Matrix_to_vector(many_S[,,ord[n]]),
                      Matrix_to_vector(many_S[,,ord[n-1]]),
                      Matrix_to_vector(many_S[,,ord[n-2]]),
                      Matrix_to_vector(many_S[,,ord[n-3]]),
                      Matrix_to_vector(many_S[,,ord[n-4]]))


short_sol = (sweep(suggestedSol,2,colMeans(encoded_train),"-"))%*%(loading)

### Number of Principal Components
k = 7



lower_1 = short_sol[1,1:k] - 0.5
upper_1 = short_sol[1,1:k] + 0.5

GA1 <- ga(type = "real-valued", 
          fitness =  function(x) test_point_depth(x,D, many_S),
          lower = lower_1, 
          upper = upper_1, 
          suggestions = short_sol[,1:k],
          popSize = 50, maxiter = 30)




ga_sol_unreduced = as.vector(loading[,1:length(GA1@solution)]%*%as.vector(GA1@solution) + colMeans(encoded_train))


vector_to_matrix = function(Cholvec){
  dim = sqrt(2*length(Cholvec)+(1/4)) - 0.5
  z = matrix(0,dim,dim)
  z[upper.tri(z, diag = TRUE)] <- Cholvec
  return(t(z)%*%z)
}

GA_So = vector_to_matrix(ga_sol_unreduced)


diagS <- diag(GA_So)
GA_Corr = diag(1/sqrt(diagS))%*%GA_So%*%diag(1/sqrt(diagS))


Diff <- array(0, dim = c(p, p, 2))


Diff[, , 1] <- GA_Corr
Diff[, , 2] <- diag(p)

out_of_sample_error = CovDist(Diff, method = "AIRM")[1,2]
round(out_of_sample_error, 3)

Diff <- array(0, dim = c(p, p, 2))


Diff[, , 1] <- many_S[,,ord[n]]
Diff[, , 2] <- diag(p)

in_sample_error = CovDist(Diff, method = "AIRM")[1,2]
round(in_sample_error, 3)


