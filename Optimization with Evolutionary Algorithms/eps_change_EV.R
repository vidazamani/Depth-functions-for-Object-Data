# changing eps to witness the impact of adding outliers 

library(ICtest)
library(abind)
library(CovTools)
library(ggplot2)
library(dplyr)
library(GA)
library(shapes)


##

set.seed(2222)


# Tukey's depth for object data (Halfspace metric depth)
MHD <- function(D){
  
  n <- nrow(D)
  p = matrix(0,n,n)
  tu = c()
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      s = 0
      for (k in 1:n) {
        if (D[k,i]<= D[k,j]){
          s=s+1
        }
        p[i,j] = (1/n)*s
        
      }
    }
  }
  
  for (y in 1:n) {
    Q = NULL
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        
        if (D[y,i]<=D[y,j]){
          Q = c(Q,p[i,j])
        }
      }
    }
    tu[y] = min(Q)
  }
  
  
  
  return(tu)
  
}



# Metric Spatial depth

MSD <- function(D){
  
  n <- nrow(D)
  
  res <- rep(0, n)
  
  for(k in 1:n){
    
    res_now <- 0
    
    d_row <- D[k, ]
    
    index_set <- setdiff(1:n, k)
    
    index_set <- setdiff(index_set, which(d_row < 1e-6))
    
    for(i in index_set){
      
      for(j in index_set){
        
        temp <- d_row[i]/d_row[j]
        
        res_now <- res_now + temp + 1/temp - D[i, j]^2/(d_row[i]*d_row[j])
        
      }
      
    }
    
    res[k] <- res_now
    
  }
  
  res <- res/n^2
  
  1 - 0.5*res
  
}

# Lens depth for Object data

EMD <- function(D){
  
  n <- nrow(D)
  maximum = c()
  lens_depth = c()
  
  for (p in 1:n) {
    s = 0
    for (i in 1:n) {
      for (j in 1:n) {
        
        if (i<j) {
          maximum[i] = max(D[i,p],D[j,p])
          if (D[i,j]>maximum[i]) {
            s = s+1
          }
          
        }
        
      }
      
    }
    
    
    lens_depth[p]  = (1/choose(n,2))*s
    
  }
  
  return(lens_depth)
  
}


# Oja depth second version for object data

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
    p[p<0] = 0 ### max(0,B) in paper
    area[k] = sum(sqrt(p))
    ojadepth[k] = 1 / (1 + (1/(0.5*(n^2-n))*(area[k])))
  }
  
  return(ojadepth)
  
}

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

vector_to_matrix = function(Cholvec){
  dim = sqrt(2*length(Cholvec)+(1/4)) - 0.5
  z = matrix(0,dim,dim)
  z[upper.tri(z, diag = TRUE)] <- Cholvec
  return(t(z)%*%z)
}

############## Change eps
p <- 5

# A total of n = 100 random covariance matrices
n <- 100

# Proportion of outliers
eps <- c(0.05)

# Non-outlying sample size
n0 <- floor((1 - eps)*n)

# Outlying sample size
n1 <- n - n0

# How outlying the outliers are
mu = 3

# Iteration 
iter = 10


D = array(0, dim = c(n, n, iter))
tudeepest = array(0, dim = c(p, p, iter))
spdeepest = array(0, dim = c(p, p, iter))
lensdeepest = array(0, dim = c(p, p, iter))
Oja_GA_deepest = array(0, dim = c(p, p, iter))
GA_So = array(0, dim = c(p, p, iter))
ojav2deepest = array(0, dim = c(p, p, iter))

MHDD = c(); MSDD = c() ; EMDD = c()  ; Ojav2D = c() ; diagSS = c()



tudiff1 = c() ; spdiff1 = c() ; lendiff1 = c() ; out_of_sample_error = c() ; ojav2diff1 = c()
tumean = c(); spmean = c() ; lenmean = c() ; ojaEV = c() ; ojav2mean = c()



for (t in 1:2) {
  # Proportion of outliers
  eps[t+1] <- eps[t] + 0.02
  
  # Non-outlying sample size
  n0[t] <- floor((1 - eps[t])*n)
  
  # Outlying sample size
  n1[t] <- n - n0[t]
  
  for (j in 1:iter) {  
    many_S_non_outliers <- replicate(n0[t], {
      
      U <- rorth(p)
      D <- diag(exp(rnorm(p)))
      S <- U%*%D%*%t(U)
      
      diagS <- diag(S)
      diag(1/sqrt(diagS))%*%S%*%diag(1/sqrt(diagS))
      
    })
    
    
    many_S_outliers <- replicate(n1[t], {
      
      U <- rorth(p)
      D <- diag(exp(rnorm(p, mu)))
      S <- U%*%D%*%t(U)
      
      diagS <- diag(S)
      diag(1/sqrt(diagS))%*%S%*%diag(1/sqrt(diagS))
      
    })
    
    many_S <- abind(many_S_non_outliers, many_S_outliers)
    
    D[,,j] = CovDist(many_S, method = "AIRM")
    
    MHDD[j] = order(MHD(D[,,j]))[n] ; MSDD[j] = order(MSD(D[,,j]))[n] 
    EMDD[j] = order(EMD(D[,,j]))[n] ; Ojav2D[j] = order(Oja_v2(D[,,j]))[n]
    
    
    
    tudeepest[,,j] = many_S[,,MHDD[j]]
    spdeepest[,,j] = many_S[,,MSDD[j]]
    lensdeepest[,,j] = many_S[,,EMDD[j]]
    ojav2deepest[,,j] = many_S[,,Ojav2D[j]]
    
    
    ### GA on Oja depth
    lower_1 = Matrix_to_vector(ojav2deepest[,,j]) - 0.5
    upper_1 = Matrix_to_vector(ojav2deepest[,,j]) + 0.5
    
    GA <- ga(type = "real-valued", 
              fitness =  function(x) test_point_depth(x,D[,,j], many_S),
              lower = lower_1, 
              upper = upper_1, 
              suggestions = Matrix_to_vector(ojav2deepest[,,j]),
              popSize = 10, maxiter = 10)
    

    
    GA_So[,,j] = vector_to_matrix(GA@solution)
    diagSS <- diag(GA_So[,,j])
    Oja_GA_deepest[,,j] = diag(1/sqrt(diagSS))%*%GA_So[,,j]%*%diag(1/sqrt(diagSS))
    
    
    
    Diff <- array(0, dim = c(p, p, 2))
    
    Diff[, , 1] <- Oja_GA_deepest[,,j]
    Diff[, , 2] <- diag(p)
    
    out_of_sample_error[j] = CovDist(Diff, method = "AIRM")[1,2]
    
    
    Diff[, , 1] <- tudeepest[,,j]
    Diff[, , 2] <- diag(p)
    
    tudiff1[j] = CovDist(Diff, method = "AIRM")[1,2]
    
    
    
    Diff[, , 1] <- spdeepest[,,j]
    spdiff1[j] = CovDist(Diff, method = "AIRM")[1,2]
    
    
    
    Diff[, , 1] <- lensdeepest[,,j]
    lendiff1[j] = CovDist(Diff, method = "AIRM")[1,2]
    
    
    Diff[, , 1] <- ojav2deepest[,,j]
    ojav2diff1[j] = CovDist(Diff, method = "AIRM")[1,2]
    
    
  }
  
  tumean[t] = mean(tudiff1) ; spmean[t] = mean(spdiff1) ; lenmean[t] = mean(lendiff1)
  ojaEV[t] = mean(out_of_sample_error)  ; ojav2mean[t] = mean(ojav2diff1)  
}





data_eps <- data.frame(x = eps[-(t+1)], a = tumean, b = spmean,
                       q = lenmean , h = ojaEV , z = ojav2mean)




ggplot(data_eps, aes(x)) +
  geom_line(aes(y = a, color = "Tukey"), lwd = 1) +
  geom_line(aes(y = b, color = "Spatial"), lwd = 1) +
  geom_line(aes(y = q, color = "Lens"), lwd = 1)+
  geom_line(aes(y = h, color = "Oja depth optimized with EV"), lwd = 1)+
  geom_line(aes(y = z, color = "Oja depth 2nd version"), lwd = 1)+
  scale_color_manual(values = c("Tukey" = "blue", "Spatial" = "coral2", "Lens" = "green",
                                "Oja depth optimized with EV" = "purple" ,
                                "Oja depth 2nd version" = "cyan3")) +
  guides(color = guide_legend(title = "Depth Method")) +
  ggtitle("The impact of adding outliers") +
  xlab("eps (Proportion of outliers)") +
  ylab("mean of error")


round(data_eps %>% summarize(mean(a), mean(b), mean(q), mean(h) , mean(z)),4)






