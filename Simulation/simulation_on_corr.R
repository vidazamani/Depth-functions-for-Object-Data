library(ICtest)
library(abind)
library(CovTools)
library(ggplot2)
library(dplyr)
library(purrr)
library(parallel)
library(dplyr)
library(viridis)
library(ggsci)


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



source("C:/Users/vizama/Documents/1st paper/Box_Parallel/Codes/All Metric depth functions.R")


simulation_with_corr <- function(p,n,eps,mu){
  
  
  # generate data
  many_S <- gener_corr_data(p,n,eps,mu)
  
  # Compute the distance matrix between the corr matrices
  
  
  D <- CovDist(many_S, method = "AIRM")
  
  
  ############## Each depth running time ###################
  
  start_time = Sys.time()
  DMHD <- MHD(D)
  end_time = Sys.time()
  MHDt <- end_time - start_time
  
  start_time = Sys.time()
  DMSD <- MSD(D)
  end_time = Sys.time()
  MSDt <- end_time - start_time
  
  
  start_time = Sys.time()
  DMLD <- MLD(D)
  end_time = Sys.time()
  MLDt <- end_time - start_time
  
  
  start_time = Sys.time()
  MOD2D <- MOD2(D)
  end_time = Sys.time()
  MOD2t <- end_time - start_time
  
  
  start_time = Sys.time()
  MOD3D <- MOD3(D)
  end_time = Sys.time()
  MOD3t <- end_time - start_time
  
  ###########################################
  
  ##### Computing Error #####
  
  MHDD = order(DMHD)[n] ; MSDD = order(DMSD)[n] 
  MLDD = order(DMLD)[n] ; MODD2 = order(MOD2D)[n]
  MODD3 = order(MOD3D)[n]
  
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
  
  
  depths_error <- c(tudiff1, spdiff1, lendiff1, oja2diff1, oja3diff1)
  time_consumption <- c(MHDt,MSDt,MLDt,MOD2t,MOD3t)
  
  result <- data.frame(depths_error,time_consumption)
  return(result)
  
}

starttime <- Sys.time()


######### Parameters 

# number of random corr matrices
n <- c(10,20,30,40,50,60)

# Matrix Dimension
p = c(3, 5, 10)

# Proportion of outliers (To check the robustness of our methods)
eps = c(0.05, 0.30)


# How outlying the outliers are
mu = 3



par_grid <- expand.grid(sample_size = n,
                        matrix_dimension = p,
                        outlier_rate = eps, 
                        mean = mu)


par_n <- nrow(par_grid)



iter = 200
temp <- list()
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


for (j in 1:iter) {
  
  set.seed(seed_vec[j])
  
  temp <- parSapply(cl, 1:par_n, function(i) simulation_with_corr(par_grid$matrix_dimension[i],
                                                                  par_grid$sample_size[i],
                                                                  par_grid$outlier_rate[i],
                                                                  par_grid$mean[i]))
  
  output <- matrix(unlist(temp), ncol = 5, byrow = TRUE)
  
  a = lapply(1:5, function(i) data.frame(par_grid, 
                                         errorD = output[seq(1,2*length(n)*length(p)*length(eps),2),i],
                                         howlong = output[seq(2,2*length(n)*length(p)*length(eps),2),i],
                                         method = i))
  
  b = list_rbind(a)
  b_all <- rbind(b_all,b)
  
}





b_all$method <- factor(b_all$method, 1:5, c("Metric Half-space depth",
                                            "Metric Spatial depth",
                                            "Metric Lens depth",
                                            "Metric Oja depth 2D",
                                            "Metric Oja depth 3D"))

bagg = b_all %>% 
  group_by(method,mean,outlier_rate,matrix_dimension,sample_size) %>% 
  summarise(avg_error = mean(errorD), avg_time = mean(howlong))

dput(bagg,file = "C:\\Users\\vizama\\Documents\\1st paper\\Box_Parallel\\data\\corr_sim_data.txt")



stopCluster(cl)


endtime <- Sys.time()
duration <- endtime - starttime



############# Visualization ###############

corr_sim_data <- dget("C:\\Users\\vizama\\Documents\\1st paper\\Box_Parallel\\data\\corr_sim_data.txt")

### 2 * 3 panel

pdf(file = 
    "C:/Users/vizama/Documents/1st paper/Box_Parallel/pic/Avg_Error/original.pdf",  
    width = 10, # The width of the plot in inches
    height = 4) # The height of the plot in inches


# New facet label names for matrix dimension variable
mtd.labs <- c("p = 3", "p = 5", "p = 10")
names(mtd.labs) <- c(3, 5, 10)

# New facet label names for outlier rate variable
otl.labs <- c("Outlier rate = 5%", "Outlier rate = 30%")
names(otl.labs) <- c("0.05", "0.3")



plot = ggplot(corr_sim_data,aes(x = sample_size, y = avg_error, col = method))+ 
  facet_grid(outlier_rate~matrix_dimension,
             labeller = labeller(matrix_dimension = mtd.labs,
                                 outlier_rate = otl.labs))+
  geom_line(linewidth = 0.7)+
  scale_x_continuous(name="Sample size")+
  scale_y_log10(name = 'Estimation error')+
  theme_bw()+ 
  theme(plot.title = element_text(size = 12), 
        plot.subtitle = element_text(size = 9),
        legend.position="bottom",
        legend.title = element_blank())





plot+scale_color_brewer(palette = 'Accent', name = 'Method')

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot+scale_color_manual(values = cbp1)

plot+scale_color_brewer(palette = 'Paired')

plot+scale_color_brewer(palette = 'Set1')

dev.off()





#### plot with title and subtitle 
# plot+labs(title = "Performance of each Metric Depth Function", 
#        subtitle = "when matrix dimnesion (P) and outlier rate increase")



#### 3 * 2 panel


pdf(file = "C:/1stPaper/Codes/Results_first_simulation/Pics/Avg_Error/blind32.pdf",  
    width = 8, # The width of the plot in inches
    height = 5) # The height of the plot in inches


# New facet label names for matrix dimension variable
mtd.labs <- c("p = 3", "p = 5", "p = 10")
names(mtd.labs) <- c(3, 5, 10)

# New facet label names for outlier rate variable
otl.labs <- c("Outlier rate = 5%", "Outlier rate = 30%")
names(otl.labs) <- c("0.05", "0.3")



plot = ggplot(bagg,aes(x = sample_size, y = avg_error, col = method))+ 
  facet_grid(matrix_dimension~outlier_rate,
             labeller = labeller(matrix_dimension = mtd.labs,
                                 outlier_rate = otl.labs))+
  geom_line(linewidth = 0.7)+
  scale_x_continuous(name="Sample size")+
  scale_y_log10(name = 'Estimation error')+
  theme_bw()+ 
  theme(plot.title = element_text(size = 12), 
        plot.subtitle = element_text(size = 9),
        legend.position="bottom",
        legend.title = element_blank())





plot+scale_color_brewer(palette = 'Accent', name = 'Method')

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot+scale_color_manual(values = cbp1)

plot+scale_color_brewer(palette = 'Paired')

plot+scale_color_brewer(palette = 'Set1')

dev.off()






##### Running time plot

### 2 * 3 panel

pdf(file = 
      "C:/Users/vizama/Documents/1st paper/Box_Parallel/pic/Avg_Time/accent T.pdf",  
    width = 10, # The width of the plot in inches
    height = 4) # The height of the plot in inches


# New facet label names for matrix dimension variable
mtd.labs <- c("p = 3", "p = 5", "p = 10")
names(mtd.labs) <- c(3, 5, 10)

# New facet label names for outlier rate variable
otl.labs <- c("Outlier rate = 5%", "Outlier rate = 30%")
names(otl.labs) <- c("0.05", "0.3")


plot2 = ggplot(corr_sim_data,aes(x = sample_size, y = avg_time, col = method))+ 
  facet_grid(outlier_rate~matrix_dimension,
             labeller = labeller(matrix_dimension = mtd.labs,
                                 outlier_rate = otl.labs))+
  geom_line(linewidth = 0.7)+
  scale_y_log10(name = 'log of running time')+
  scale_x_continuous(name="sample size")+
  theme_bw()+ 
  theme(plot.title = element_text(size = 12), 
        plot.subtitle = element_text(size = 9),
        legend.position="bottom",
        legend.title = element_blank())


plot2+scale_color_brewer(palette = 'Accent', name = 'Method')

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot+scale_color_manual(values = cbp1)

plot+scale_color_brewer(palette = 'Paired')

plot+scale_color_brewer(palette = 'Set1')

dev.off()

#####

plot3 = ggplot(corr_sim_data,aes(x = sample_size, y = avg_error, col = method,linetype = factor(outlier_rate)))+ 
  facet_grid(matrix_dimension~.)+
  geom_line(linewidth = 0.7)+
  scale_x_continuous(name="sample size")+
  theme_bw()+ 
  labs(title = "Performance of each Metric Depth Function", 
       subtitle = "when matrix dimnesion is 3 and 5% of Distribution is contaminated",
       y = 'Bias in log scale')+ 
  theme(plot.title = element_text(size = 12), 
        plot.subtitle = element_text(size = 9))+
  scale_y_log10()


###### Uni Plot


plotuni = ggplot(p3eps0.3,aes(x = sample_size, y = avg_error, col = method))+ 
  geom_line(linewidth = 0.7)+
  scale_x_continuous(name="sample size")+
  scale_y_continuous(name = 'Average Error')+
  theme_bw()+ 
  labs(title = "Performance of each Metric Depth Function", 
       subtitle = "when matrix dimnesion is 3 and 5% of Distribution is contaminated")+ 
  theme(plot.title = element_text(size = 12), 
        plot.subtitle = element_text(size = 9))















