library(ICtest)
library(shapes)
library(factoextra)
library(CovTools)
library(GA)
library(abind)
library(ggplot2)
library(dplyr)
library(purrr)
library(parallel)
library(dplyr)
library(viridis)
library(ggsci)
library(ggpubr)


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
    D <- diag(exp(c(rnorm(1, mu), rnorm(p - 1, -mu))))
    S <- U%*%D%*%t(U)
    
    diagS <- diag(S)
    diag(1/sqrt(diagS))%*%S%*%diag(1/sqrt(diagS))
    
    
  })
  
  many_S <- abind(many_S_non_outliers, many_S_outliers)
  return(many_S)
}



# Oja-2D depth between train sample and test point
MOD2TRTS <- function(D, d){
  
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
  p[p< 1e-6] = 0
  area = sum(sqrt(p))
  ojadepth = 1 / (1 + (1/(0.5*(n^2-n))*(area)))
  
  return(ojadepth)
}


# Helper functions

# Converts a correlation matrix into its encoding
matrix_to_vector = function(CovMatrix){
  U <- chol(CovMatrix)
  U[upper.tri(U, diag = TRUE)]
}


# Converts an encoding into a correlation matrix
vector_to_matrix = function(Cholvec){
  dim = sqrt(2*length(Cholvec)+(1/4)) - 0.5
  z = matrix(0,dim,dim)
  z[upper.tri(z, diag = TRUE)] <- Cholvec
  return(t(z)%*%z)
}


source("C:/Users/vizama/Documents/1st paper/Box_Parallel/Codes/All Metric depth functions.R")



GA_oja2 = function(n, p, eps, mu, iterga, PEL){
  
  
  # generate data
  many_S <- gener_corr_data(p,n,eps,mu)
  
  start_timepca = Sys.time()
  
  # Compute the distance matrix between the corr matrices
  D <- CovDist(many_S, method = "AIRM")
  
  # Corresponding encoding length
  p_enc <- 0.5*p*(p + 1)
  
  
  
  #### Some auxiliary objects needed by PCA
  
  # Encoded training data
  encoded_train <- matrix(0, n, p_enc)
  
  for (i in 1:n) {
    encoded_train[i, ] <- matrix_to_vector(many_S[, , i])
  }
  
  
  ###### NEW LINE ####################
  # Whether we want to center the principal components or not
  centering <- FALSE
  
  # Mean vector of the encoded training data
  
  if(centering){
    mean_vec <- colMeans(encoded_train)
  }
  
  if(!centering){
    mean_vec <- rep(0, p_enc)
  }
  
  
  
  # PCA either with centering or without
  res.pca <- prcomp(encoded_train, center = centering, scale = FALSE)
  
  
  # Loadings
  loading <- res.pca$rotation
  
  # Number of PCs is taken here as a fixed proportion of the encoding length
  num_pc <- floor(PEL*min(p_enc,n)) ###k
  
  #### PCA Fitness function
  # v = the encoded test point in the PCA-space
  # D = training data distance matrix
  # many_S = training data as the original correlation matrices
  test_point_depth = function(v, D, many_S){
    
    v_star = as.vector(loading[, 1:num_pc, drop = FALSE]%*%v + mean_vec)
    
    CM_test = vector_to_matrix(v_star) 
    
    d = c()
    
    for (i in 1:dim(many_S)[3]) {
      d[i] = distcov(many_S[, , i], CM_test)
    }
    
    MOD2TRTS(D, d)
  }
  
  # Top 5 in-sample points
  ord = order(MOD2(D), decreasing = FALSE)
  
  # Matrix of the encoded vectors of the top 5 in-sample solutions
  suggestedSol <- rbind(matrix_to_vector(many_S[, , ord[n]]),
                        matrix_to_vector(many_S[, , ord[n - 1]]),
                        matrix_to_vector(many_S[, , ord[n - 2]]),
                        matrix_to_vector(many_S[, , ord[n - 3]]),
                        matrix_to_vector(many_S[, , ord[n - 4]]))
  
  # The above converted to the principal components
  
  short_sol <- (sweep(suggestedSol, 2, mean_vec, "-"))%*%(loading[, 1:num_pc, drop = FALSE])
  
  # The search interval
  # (HOW SHOULD THIS BE CHOSEN? IS +- 0.5 GOOD?)
  lower_1 = short_sol[1, 1:num_pc] - 0.5
  upper_1 = short_sol[1, 1:num_pc] + 0.5
  
  
  # The actual genetic algorithm for optimization
  GA1 <- ga(type = "real-valued", 
            fitness =  function(x) test_point_depth(x, D, many_S),
            lower = lower_1, 
            upper = upper_1, 
            suggestions = short_sol,
            popSize = 500, maxiter = iterga)
  
  
  # The found solution in the original space
  ga_sol_unreduced <- as.vector(loading[, 1:num_pc, drop = FALSE]%*%as.vector(GA1@solution) + mean_vec)
  
  # The found solution as a matrix
  GA_So = vector_to_matrix(ga_sol_unreduced)
  
  # The result is not necessary a correlation matrix since the 
  # optimization might cause us to leave the space of encoded correlation matrices and enter the space
  # of encoded covariance matrices 
  # Thus the result needs to be converted into a correlation matrix
  diagS <- diag(GA_So)
  GA_Corr = diag(1/sqrt(diagS))%*%GA_So%*%diag(1/sqrt(diagS))
  
  # Find the error between the solution found by GA and the true deepest point
  Diff <- array(0, dim = c(p, p, 2))
  
  Diff[, , 1] <- GA_Corr ## best solution based on GA
  Diff[, , 2] <- diag(p)
  
  out_of_sample_error = CovDist(Diff, method = "AIRM")[1,2]
  
  End_timepca = Sys.time()
  
  pcatime = difftime(End_timepca ,start_timepca, units = 'secs')
  
  # Compare this to the error between the deepest in-sample point and the true deepest point
  
  # set.seed(1234)
  
  Diff <- array(0, dim = c(p, p, 2))
  
  # generate data
  many_So <- gener_corr_data(p,n,eps,mu)
  
  start_timeo = Sys.time()
  
  # Compute the distance matrix between the corr matrices
  D <- CovDist(many_So, method = "AIRM")
  ordd = order(MOD2(D), decreasing = FALSE)
  
  Diff[, , 1] <- many_So[,,ordd[n]]
  Diff[, , 2] <- diag(p)
  
  in_sample_error = CovDist(Diff, method = "AIRM")[1,2]
  
  end_timeo = Sys.time()
  
  Oja2time = difftime(end_timeo, start_timeo, units = 'secs')
  
  return(data.frame(OSE = out_of_sample_error,
                    ISE = in_sample_error,pcatime,Oja2time))
  
}

######### Parameters 

# Matrix Dimension
p <- 5


# number of random corr matrices
n <- c(10,20,30,40,50,60)
# n = 10

# How outlying the outliers are
mu = 3

# Proportion of outliers (To check the robustness of our methods)
eps = 0.3


iterga = 30

# proportion of the encoding length
PEL = c(0.25, 0.5, 0.75)
# PEL = 0.25

par_grid <- expand.grid(sample_size = n,
                        matrix_dimension = p,
                        outlier_rate = eps, 
                        mean = mu,
                        GA_iteration = iterga,
                        PCA_n = PEL)


par_n <- nrow(par_grid)

repli = 200
temp <- list()
output <- matrix(0)
a <- list()
b <- data.frame()
b_all <- data.frame()

begin <- Sys.time()

cl <- makeCluster(4)
clusterExport(cl, ls(), envir = environment())
clusterExport(cl, c("rorth","abind", "CovDist", 'ga', 'distcov'),
              envir = environment())


clusterEvalQ(cl, set.seed(2222))
set.seed(1111)
seed_vec <- sample(1:100000, repli)

for (j in 1:repli) {
  
  set.seed(seed_vec[j])
  temp <- parSapply(cl, 1:par_n,
                    function(i) GA_oja2(par_grid$sample_size[i],
                                        par_grid$matrix_dimension[i],
                                        par_grid$outlier_rate[i],
                                        par_grid$mean[i],
                                        par_grid$GA_iteration[i],
                                        par_grid$PCA_n[i]))
  
  output <- matrix(unlist(temp), ncol = 1, byrow = TRUE)
  
  a = lapply(1:1, function(i) data.frame(par_grid, 
                                         OSE = output[seq(1,4*length(n)*length(PEL),4),i],
                                         ISE = output[seq(2,4*length(n)*length(PEL),4),i],
                                         pcatime = output[seq(3,4*length(n)*length(PEL),4),i],
                                         oja2time = output[seq(4,4*length(n)*length(PEL),4),i]))
  
  b = list_rbind(a)
  b_all <- rbind(b_all,b)
  
}

ending <- Sys.time()
runningtime <- ending - begin

insdata = b_all
insdata$PCA_n <- as.character(insdata$PCA_n)

outdata = b_all
outdata$PCA_n <- as.character(outdata$PCA_n)

timedata = b_all
timedata$PCA_n <- as.character(timedata$PCA_n)


data1 <- insdata %>% 
  group_by(PCA_n,sample_size) %>% 
  summarise(avg_in_sample_error = mean(ISE),
            avg_outof_sample_error = mean(OSE))

data1 = data1[1:length(n),]
data1 <-  data1[c(2,3)]



data2 = outdata %>% 
  group_by(PCA_n,sample_size) %>% 
  summarise(avg_in_sample_error = mean(ISE),
            avg_outof_sample_error = mean(OSE))

data3 <- timedata %>% 
  group_by(PCA_n,sample_size) %>% 
  summarise(avg_pca_time = mean(pcatime),
            avg_oja_time = mean(oja2time))


data4 <-  data3[c(2,4)] %>% group_by(sample_size) %>% 
  summarise(avg_oja_time = mean(avg_oja_time))


###### Estimation error plot

p1 <- ggplot() +  
  geom_line(data = data2, aes(x = sample_size,
                              y = avg_outof_sample_error,
                              linetype = PCA_n,
                              colour = PCA_n),
            linewidth = 1) +  
  geom_line(data = data1, aes(x = sample_size,
                              y = avg_in_sample_error,
                              linetype = 'In-Sample Estimator',
                              colour = 'In-Sample Estimator'),
            linewidth = 1)+
  scale_linetype_manual(values = c("0.25"= "dashed",
                                   "0.5" = "solid",
                                   "0.75" = "dotted",
                                   "In-Sample Estimator" = "dotdash"))+
  scale_color_manual(values = c("0.25"= "blue",
                                "0.5" = "green",
                                "0.75" = "purple",
                                "In-Sample Estimator" = "red"))+
  scale_x_continuous(name="Sample Size")+
  scale_y_continuous(name = 'Estimation Error')+
  theme_bw()+ 
  theme(plot.title = element_text(size = 12), 
        plot.subtitle = element_text(size = 9),
        legend.position="bottom",
        legend.title = element_blank())



#### Time plot in logaritmic scale

p2 <- ggplot() +  
  geom_line(data = data3, aes(x = sample_size,
                              y = avg_pca_time,
                              linetype = PCA_n,
                              colour = PCA_n),
            linewidth = 1) +  
  geom_line(data = data4, aes(x = sample_size,
                              y = avg_oja_time,
                              linetype = 'In-Sample Estimator',
                              colour = 'In-Sample Estimator'),
            linewidth = 1)+
  scale_linetype_manual(values = c("0.25"= "dashed",
                                   "0.5" = "solid",
                                   "0.75" = "dotted",
                                   "In-Sample Estimator" = "dotdash"))+
  scale_color_manual(values = c("0.25"= "blue",
                                "0.5" = "green",
                                "0.75" = "purple",
                                "In-Sample Estimator" = "red")) +
  scale_x_continuous(name="Sample Size")+
  scale_y_log10(name = 'Estimation time in sec')+
  theme_bw()+ 
  theme(plot.title = element_text(size = 12), 
        plot.subtitle = element_text(size = 9),
        legend.position="bottom",
        legend.title = element_blank())+
  guides(color = guide_legend(title = "PCA prop"),
         linetype = guide_legend(title = "PCA prop"))






#### Time plot 

p3 <- ggplot() +  
  geom_line(data = data3, aes(x = sample_size,
                              y = avg_pca_time,
                              linetype = PCA_n,
                              colour = PCA_n),
            linewidth = 1) +  
  geom_line(data = data4, aes(x = sample_size,
                              y = avg_oja_time,
                              linetype = 'In-Sample Estimator',
                              colour = 'In-Sample Estimator'),
            linewidth = 1)+
  scale_linetype_manual(values = c("0.25"= "dashed",
                                   "0.5" = "solid",
                                   "0.75" = "dotted",
                                   "In-Sample Estimator" = "dotdash"))+
  scale_color_manual(values = c("0.25"= "blue",
                                "0.5" = "green",
                                "0.75" = "purple",
                                "In-Sample Estimator" = "red")) +
  scale_x_continuous(name="Sample Size")+
  scale_y_continuous(name = 'Estimation time')+
  theme_bw()+ 
  theme(plot.title = element_text(size = 12), 
        plot.subtitle = element_text(size = 9),
        legend.position="bottom",
        legend.title = element_blank())



stopCluster(cl)

pdf("C:/1stPaper/Codes/GA/ga plot.pdf", 6 ,4)
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()

