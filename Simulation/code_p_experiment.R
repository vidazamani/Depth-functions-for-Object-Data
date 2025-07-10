# Oja-paper distance computation simulation


library(tidyverse)
library(MetricDepth)



rsphere <- function(n, p, lambda){
  x <- matrix(rnorm(n*p, lambda, 1), n, p)
  x_norm <- apply(x, 1, function(v) sqrt(sum(v^2)))
  return(sweep(x, 1, x_norm, "/"))
}



cont_sample <- function(n, p, eps, lambda){
  n1 <- floor((1 - eps)*n)
  n2 <- n - n1
  
  x1 <- rsphere(n1, p, lambda)
  x2 <- rsphere(n2, p, -1)
  
  rbind(x1, x2)
}



single_simu <- function(n, p){
  x <- cont_sample(n, p, 0.10, 5)
  
  suppressWarnings(D <- acos(tcrossprod(x)))
  diag(D) <- rep(0, nrow(D))
  
  DMHD <- MHD_cpp(D)
  DMSD <- MSD_cpp(D)
  DMLD <- MLD_cpp(D)
  MOD2D <- MOD2_cpp(D)
  MOD3D <- MOD3_cpp(D)
  
  MHDD <- order(DMHD)[n]
  MSDD <- order(DMSD)[n] 
  MLDD <- order(DMLD)[n]
  MODD2 <- order(MOD2D)[n]
  MODD3 <- order(MOD3D)[n]
  
  RAND <- x[sample(n, 1), ]

  
  tudiff1 = acos(sum(x[MHDD, ]*(1/sqrt(p)*rep(1,p))))
  spdiff1 = acos(sum(x[MSDD, ]*(1/sqrt(p)*rep(1,p))))
  lendiff1 = acos(sum(x[MLDD, ]*(1/sqrt(p)*rep(1,p))))
  oja2diff1 = acos(sum(x[MODD2, ]*(1/sqrt(p)*rep(1,p))))
  oja3diff1 = acos(sum(x[MODD3, ]*(1/sqrt(p)*rep(1,p))))
  randdiff1 <- acos(sum(RAND*(1/sqrt(p)*rep(1,p))))
  
  # Experimental part
  # z <- colMeans(x[which(DMHD == max(DMHD)), , drop = FALSE])
  # z <- z/sqrt(sum(z^2))
  # tudiff1 = acos(sum(z*(1/sqrt(p)*rep(1,p))))
  
  
  data.frame(n = n, p = p, method = c(1:6), res = c(tudiff1, spdiff1, lendiff1, oja2diff1, oja3diff1, randdiff1))
}


# n <- 100
# p <- 5
# 
# single_simu(n, p)



# Batches 1-5: n = 100
# Batch 6: n = 150
# Batch 7: n = 200
# Always 20 reps per batch

# Next up: batch 11
# Batch 11 run with set.seed(30111988 + batch - 1 + 100)

n <- 200
p_set <- c(50, 100, 200, 400, 800, 1600)

reps_per_set <- 20

batch <- 11

set.seed(30111988 + batch - 1)
res <- NULL

for(p in p_set){
  for(j in 1:reps_per_set){
    temp <- single_simu(n, p)
    res <- bind_rows(res, temp)
    if(j %% 10 == 0){
      print(j)
    }
  }
  print(p)
}

  
  
write.table(res, paste0("/Users/jomivi/Library/Mobile Documents/com~apple~CloudDocs/Work/Supervision/Vida/Paper_1_revision/results_p_experiment/res_", batch, ".txt"))


# res

res <- NULL

for(i in c(1:5, 7:11)){
  # temp <- read.table(paste0("C:\\Users\\joniv\\iCloudDrive\\Work\\Supervision\\Lauri\\Paper_2\\Revision\\Timing_results\\res_", i, ".txt"))
  temp <- read.table(paste0("/Users/jomivi/Library/Mobile Documents/com~apple~CloudDocs/Work/Supervision/Vida/Paper_1_revision/results_p_experiment/res_", i, ".txt"))
  res <- bind_rows(res, temp)
}




res_agg <- res %>%
  group_by(method, n, p) %>%
  summarize(avg_res = mean(res)) %>%
  filter(method != 6)

res_rand_1 <- structure(list(method = c(6, 6, 6, 6, 6, 6),
                             p = c(50, 100, 200, 400, 800, 1600),
                             n = c(100, 100, 100, 100, 100, 100),
                             avg_res = c(0.408273612481363, 0.416033071123464, 0.402427058769401, 0.417530539066569, 0.420584191113979, 0.410348393627898)),
                        class = "data.frame", row.names = c(NA, -6L))

res_rand_2 <- structure(list(method = c(6, 6, 6, 6, 6, 6), p = c(50, 100, 200, 
                                                                 400, 800, 1600), n = c(200, 200, 200, 200, 200, 200), avg_res = c(0.410778443867579, 
                                                                                                                                   0.409466705107033, 0.402345852895837, 0.404709226888634, 0.415121230864784, 
                                                                                                                                   0.402763411742536)), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                            -6L))


# res_agg <- bind_rows(res_agg, res_rand) %>%
#   mutate(method = factor(method, levels = 1:6, labels = c("MHD", "MSD", "MLD", "MOD2", "MOD3", "RAND")))

res_agg <- bind_rows(res_agg, res_rand_1, res_rand_2) %>%
  mutate(method = factor(method, levels = 1:6, labels = c("MHD", "MSD", "MLD", "MOD2", "MOD3", "RAND"))) %>%
  mutate(n = factor(n))


pdf(paste0("/Users/jomivi/Library/Mobile Documents/com~apple~CloudDocs/Work/Supervision/Vida/Paper_1_revision/plot_p_experiment.pdf"), width = 5, height = 4)

ggplot(filter(res_agg, method != "RAND"), aes(x = p, y = avg_res)) +
  geom_line(aes(col = method, linetype = n)) +
  theme_bw() +
  labs(x = "Dimension p", y = "Estimation error", col = "Method") # +
  # scale_y_continuous(transform = "log")

dev.off()

# Recompute the RAND-line

set.seed(1234)

n <- 200
p_set <- c(50, 100, 200, 400, 800, 1600)

res_rand_1 <- NULL

for(p in p_set){
  temp <- replicate(10000, {
    x <- cont_sample(n, p, 0.10, 5)
    acos(sum(x[sample(n, 1), ]*(1/sqrt(p)*rep(1,p))))
  })

  res_rand_1 <- bind_rows(res_rand_1, data.frame(method = 6, p = p, n = n, avg_res = mean(temp)))
}

dput(res_rand_1)

