# generate data function 

gener_cov_data <- function(p,n,eps,mu){
  
  # Non-outlying sample size
  n0 <- floor((1 - eps)*n)
  
  # Outlying sample size
  n1 <- n - n0
  
  many_S_non_outliers <- replicate(n0, {
    
    U <- rorth(p)
    D <- diag(exp(rnorm(p)))
    S <- U%*%D%*%t(U)
    
    S
    
  })
  
  
  many_S_outliers <- replicate(n1, {
    
    U <- rorth(p)
    D <- diag(exp(rnorm(p, mu)))
    S <- U%*%D%*%t(U)
    
    S
    
  })
  
  many_S <- abind(many_S_non_outliers, many_S_outliers)
  return(many_S)
}



source("D:/All Metric depth functions.R")

## Or
https://github.com/vidazamani/Depth-functions-for-Object-Data/blob/95d9a29cdb39ec1c33b0acc0951f93b20f8a77be/All%20Metric%20depth%20functions.R#L1





simulation_with_cov <- function(p,n,eps,mu){
  
  
  # generate data
  many_S <- gener_cov_data(p,n,eps,mu)
  
  # Compute the distance matrix between the correlation matrices
  
  
  D <- CovDist(many_S, method = "AIRM")
  
  # print("this is metric lens dapth:")
  # print(MLD(D))
  # 
  # print("this is metric half-space dapth:")
  # print(MHD(D))
  # 
  # print("this is metric spatial dapth:")
  # print(MSD(D))
  # 
  # print("this is metric Oja dapth in 2D:")
  # print(MOD2(D))
  # 
  # print("this is metric Oja dapth in 3D:")
  # print(MOD3(D))
  
  
  ##### Computing Error #####
  
  MHDD = order(MHD(D))[n] ; MSDD = order(MSD(D))[n] 
  MLDD = order(MLD(D))[n] ; MODD2 = order(MOD2(D))[n]
  MODD3 = order(MOD3(D))[n]
  
  
  tudeepest = many_S[,,MHDD]
  spdeepest = many_S[,,MSDD]
  lensdeepest = many_S[,,MLDD]
  oja2deepest = many_S[,,MODD2]
  oja3deepest = many_S[,,MODD3]
  
  
  Diff <- array(0, dim = c(p, p, 2))
  
  
  Diff[, , 1] <- tudeepest
  Diff[, , 2] <- diag(p)
  
  tudiff1 = CovDist(Diff, method = "AIRM")[1,2]
  
  
  
  Diff[, , 1] <- spdeepest
  spdiff1 = CovDist(Diff, method = "AIRM")[1,2]
  
  
  
  Diff[, , 1] <- lensdeepest
  lendiff1 = CovDist(Diff, method = "AIRM")[1,2]
  
  
  Diff[, , 1] <- oja2deepest
  oja2diff1 = CovDist(Diff, method = "AIRM")[1,2]
  
  
  Diff[, , 1] <- oja3deepest
  oja3diff1 = CovDist(Diff, method = "AIRM")[1,2]
  
  
  print(c(tudiff1, spdiff1, lendiff1, oja2diff1, oja3diff1))
  # print(tudeepest) ; print(lensdeepest)
  # print(MHDD) ; print(MSDD) ; print(MLDD)
  # print(MODD2) ; print(MODD3)
}




set.seed(2222)



######### Parameters 

p <- 5

# A total of n = 100 random covariance matrices
n <- 50

# Proportion of outliers
eps <- 0.05

# How outlying the outliers are
mu = 3




simulation_with_cov(p,n,eps,mu)

iter = 30
m = matrix(0,iter,5)
for (i in 1:iter) {
  
  m[i,] = simulation_with_cov(p,n,eps,mu)
  

}

colMeans(m)








# p_list <- c(5)
# n_list <- c(30, 40, 50, 60, 70)
# eps_list <- c(0.05)
# mu_list <- c(3)
# 
# res <- NULL
# 
# for(p in p_list){
#   for(n in n_list){
#     for(eps in eps_list){
#       for(mu in mu_list){
#         temp <- data.frame()
#         res <- rbind(res, temp)
#       }
#     }
#   }
# }


p <- 5

# A total of n = 100 random covariance matrices
n <- c(10,20)

# Proportion of outliers
eps <- 0.05

# How outlying the outliers are
mu = 3

par_grid <- expand.grid(sample_size = n, matrix_dimension = p, outlier_rate = eps, mean = mu)


par_n <- nrow(par_grid)

library(purrr)
iter = 2
b_all <- NULL
for (i in 1:iter) {
  temp <- sapply(1:par_n, function(i) simulation_with_cov(par_grid$matrix_dimension[i],
                                                          par_grid$sample_size[i],
                                                          par_grid$outlier_rate[i],
                                                          par_grid$mean[i]))
  a = lapply(1:5, function(i) data.frame(par_grid, simulationwithcov = temp[i,], method = i))
  
  b = list_rbind(a)
  b_all <- rbind(b_all,b)
}




b_all$method <- factor(b_all$method, 1:5, c("Metric Lens depth","Metric Half-space depth",
                                    "metric spatial depth",
                                    "metric Oja depth 2D",
                                    "metric Oja depth 3D"))

bagg = b_all %>% 
  group_by(method,mean,outlier_rate,matrix_dimension,sample_size) %>% summarise(avg_error = mean(simulationwithcov))


plot = ggplot(bagg,aes(x = sample_size, y = avg_error, col = method))+ geom_line()

