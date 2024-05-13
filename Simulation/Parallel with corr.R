library(ICtest)
library(abind)
library(CovTools)
library(ggplot2)
library(dplyr)
library(purrr)
library(parallel)




# generate data function 

gener_corr_data <- function(p,n,eps,mu){
  
  # Non-outlying sample size
  n0 <- floor((1 - eps)*n)
  
  # Outlying sample size
  n1 <- n - n0
  
  many_S_non_outliers <- replicate(n0, {
    
    U <- rorth(p)
    D <- diag(exp(rnorm(p)))
    S <- U%*%D%*%t(U)
    
    diagS <- diag(S)
    diag(1/sqrt(diagS))%*%S%*%diag(1/sqrt(diagS))
    
  })
  
  
  many_S_outliers <- replicate(n1, {
    
    U <- rorth(p)
    D <- diag(exp(rnorm(p, mu)))
    S <- U%*%D%*%t(U)
    
    diagS <- diag(S)
    diag(1/sqrt(diagS))%*%S%*%diag(1/sqrt(diagS))
    
    
  })
  
  many_S <- abind(many_S_non_outliers, many_S_outliers)
  return(many_S)
}



source("C:/1st Paper/Codes/All Metric depth functions.R")






simulation_with_corr <- function(p,n,eps,mu){
  
  
  # generate data
  many_S <- gener_corr_data(p,n,eps,mu)
  
  # Compute the distance matrix between the covariance matrices
  
  
  D <- CovDist(many_S, method = "AIRM")
  
  # print("this is metric lens depth:")
  # print(MLD(D))
  # 
  # print("this is metric half-space depth:")
  # print(MHD(D))
  # 
  # print("this is metric spatial depth:")
  # print(MSD(D))
  # 
  # print("this is metric Oja depth in 2D:")
  # print(MOD2(D))
  # 
  # print("this is metric Oja depth in 3D:")
  # print(MOD3(D))
  
  ############## Time consumption by each depth ###################
  
  start_time = Sys.time()
  MHD(D)
  end_time = Sys.time()
  MHDt <- end_time - start_time
  
  start_time = Sys.time()
  MSD(D)
  end_time = Sys.time()
  MSDt <- end_time - start_time
  
  
  start_time = Sys.time()
  MLD(D)
  end_time = Sys.time()
  MLDt <- end_time - start_time
  
  
  start_time = Sys.time()
  MOD2(D)
  end_time = Sys.time()
  MOD2t <- end_time - start_time
  
  
  start_time = Sys.time()
  MOD3(D)
  end_time = Sys.time()
  MOD3t <- end_time - start_time
  
  ###########################################
  
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
  
  
  # print(c(tudiff1, spdiff1, lendiff1, oja2diff1, oja3diff1))
  # print(c(MHDt,MSDt,MLDt,MOD2t,MOD3t))
  
  depths_error <- c(tudiff1, spdiff1, lendiff1, oja2diff1, oja3diff1)
  time_consumption <- c(MHDt,MSDt,MLDt,MOD2t,MOD3t)
  
  result <- data.frame(depths_error,time_consumption)
  return(result)
  # print(tudeepest) ; print(lensdeepest)
  # print(MHDD) ; print(MSDD) ; print(MLDD)
  # print(MODD2) ; print(MODD3)
}







######### Parameters 

# p <- 5
# 
# # A total of n = 100 random covariance matrices
# n <- 5
# 
# # Proportion of outliers
# eps <- 0.05
# 
# # How outlying the outliers are
# mu = 3




p <- 10

# number of random covariance matrices
n <- c(10,20,30,40,50,60)

# n <- c(10,20)

# Proportion of outliers (To check th robustness of our methods)
eps <- 0.05

# How outlying the outliers are
mu = 3

par_grid <- expand.grid(sample_size = n,
                        matrix_dimension = p,
                        outlier_rate = eps, 
                        mean = mu)


par_n <- nrow(par_grid)


# a = sapply(1:par_n, function(i) simulation_with_cov(par_grid$matrix_dimension[i],
#                                                     par_grid$sample_size[i],
#                                                     par_grid$outlier_rate[i],
#                                                     par_grid$mean[i]))
# n = c(10,20)
# 
# output <- matrix(unlist(a), ncol = 5, byrow = TRUE)
# 
# 
# 
# q = lapply(1:5, function(i) data.frame(par_grid,
#                                        errorD = output[seq(1,2*length(n),2),i],
#                                        howlong = output[seq(2,2*length(n),2),i],
#                                        method = i))




iter = 1
# temp <- matrix(0,5,par_n)
temp <- list()
# output <- matrix(0,2*length(n),5)
output <- matrix(0)
a <- list()
b <- data.frame()
b_all <- data.frame()



cl <- makeCluster(4)
clusterExport(cl, ls(), envir = environment())
clusterExport(cl, c("rorth","abind", "CovDist"), envir = environment())


clusterEvalQ(cl, set.seed(2222))
set.seed(1111)
seed_vec <- sample(1:100000, iter)

starttime <- Sys.time()

for (j in 1:iter) {
  
  set.seed(seed_vec[j])
  
  temp <- parSapply(cl, 1:par_n, function(i) simulation_with_cov(par_grid$matrix_dimension[i],
                                                                 par_grid$sample_size[i],
                                                                 par_grid$outlier_rate[i],
                                                                 par_grid$mean[i]))
  
  output <- matrix(unlist(temp), ncol = 5, byrow = TRUE)
  
  a = lapply(1:5, function(i) data.frame(par_grid, 
                                         errorD = output[seq(1,2*length(n),2),i],
                                         howlong = output[seq(2,2*length(n),2),i],
                                         method = i))
  
  b = list_rbind(a)
  b_all <- rbind(b_all,b)

}

endtime <- Sys.time()
duration <- endtime - starttime



b_all$method <- factor(b_all$method, 1:5, c("Metric Lens depth","Metric Half-space depth",
                                            "metric spatial depth",
                                            "metric Oja depth 2D",
                                            "metric Oja depth 3D"))

bagg = b_all %>% 
  group_by(method,mean,outlier_rate,matrix_dimension,sample_size) %>% 
  summarise(avg_error = mean(errorD), avg_time = mean(howlong))



plot = ggplot(bagg,aes(x = sample_size, y = avg_error, col = method))+ 
  geom_line(linewidth = 0.7)

plot2 = ggplot(bagg,aes(x = sample_size, y = avg_time, col = method))+ 
  geom_line(linewidth = 0.7)

write.table(bagg, file = "C:\\Users\\vizama\\Documents\\1st paper\\Box_Parallel\\data\\results.txt")


stopCluster(cl)
