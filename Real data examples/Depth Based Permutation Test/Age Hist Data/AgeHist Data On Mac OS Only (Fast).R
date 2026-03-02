library(HistDAWass)
library(tidyr)
library(dplyr)
library(Rcpp) ## you need to have the latest version of Rcpp
library(ggpattern) 
library(patchwork)
library(ggpubr)
library(parallel)

remotes::install_github('vidazamani/Depth-functions-for-Object-Data/MetricDepthCpp')
library(MetricDepth)


M <- attr(Age_Pyramids_2014, "M")
n_total <- nrow(M)
countries <- rownames(M)


# Create a data frame with every country and its index
df_regions <- data.frame(
  country = countries,
  index   = seq_along(countries),
  stringsAsFactors = FALSE
)


####### Europe 

western     <- c("Austria","Belgium","France","Germany",
                 "Liechtenstein","Luxembourg","Monaco",
                 "Netherlands","Switzerland", 'Ireland','United Kingdom')


southwest   <- c("Portugal","Spain",'Andorra')  


southern    <- c("Italy","Malta",'Greece',
                 "San Marino","Vatican City")


northern    <- c("Denmark","Finland","Iceland",
                 "Norway","Sweden","United Kingdom", 'Faroe Islands',
                 'Jan Mayen',
                 'Svalbard')


central     <- c("Czech Republic","Hungary","Poland",
                 "Slovakia","Slovenia")


southeastern<- c("Albania","Bosnia and Herzegovina",
                 "Bulgaria","Croatia",
                 "Montenegro","North Macedonia","Romania",
                 "Serbia","Slovenia")

eastern     <- c("Belarus","Moldova",
                 "Romania","Russia","Ukraine")

# 3. Build the two merged groups
group_west <- unique(c(western, northern, southern, southwest))
group_east <- unique(c(eastern, southeastern, central))


# Find matching indices
idx_east <- which(countries %in% group_east)
idx_west <- which(countries %in% group_west)


# Build lookup table
df_eastern <- data.frame(
  country = countries[idx_east],
  index   = idx_east,
  stringsAsFactors = FALSE
)


df_western <- data.frame(
  country = countries[idx_west],
  index   = idx_west,
  stringsAsFactors = FALSE
)




df_europe <- rbind(df_eastern,df_western)
idx_europe <- df_europe$index

############### African

african_countries <- c(
  "Algeria", "Angola", "Benin", "Botswana", "Burkina Faso", "Burundi", "Cabo Verde", "Cameroon",
  "Central African Republic", "Chad", "Comoros", "Congo (Brazzaville)", 
  "Congo (Kinshasa)", "Cote d'Ivoire", "Djibouti", "Egypt", "Equatorial Guinea", "Eritrea", 
  "Swaziland", "Ethiopia", "Gabon", "Gambia, The", "Ghana", "Guinea", "Guinea-Bissau", "Kenya",
  "Lesotho", "Liberia", "Libya", "Madagascar", "Malawi", "Mali", "Mauritania", "Mauritius",
  "Morocco", "Mozambique", "Namibia", "Niger", "Nigeria", "Rwanda",
  "Saint Helena", "Sao Tome and Principe", "Senegal", "Seychelles",
  "Sierra Leone", "Somalia", "South Africa", "South Sudan", "Sudan", "Tanzania", "Togo", 
  "Tunisia", "Uganda", "Zambia", "Zimbabwe"
)



idx_african <- which(countries %in% african_countries)


df_african <- data.frame(
  country = countries[idx_african],
  index   = idx_african,
  stringsAsFactors = FALSE
)




##### Parameters 


depth_names <- c("MOD3", "MOD2", "MHD", "MLD", "MSD")
depth_models <- list(MOD3_cpp, MOD2_cpp, MHD_cpp, MLD_cpp, MSD_cpp)




##########


full_wass_matrix <- matrix(0, n_total, n_total)

for (i in 1:(n_total - 1)) {
  
  for (j in (i + 1):n_total) {
    
    d <- WassSqDistH(
      Age_Pyramids_2014[i]@M[[1]],
      Age_Pyramids_2014[j]@M[[1]]
    )
    
    full_wass_matrix[i, j] <- sqrt(d)
    full_wass_matrix[j, i] <- sqrt(d)
  }
  
}

dim(full_wass_matrix)


compute_wass_submatrix <- function(group_idx) {
  full_wass_matrix[group_idx, group_idx, drop = FALSE]
}


pvalue_age_hist_fast <- function(iter, depth_fun) {
  
  # ---- Observed statistic ----
  
  dist_af_obs <- compute_wass_submatrix(idx_african)
  dist_eu_obs <- compute_wass_submatrix(idx_europe)
  
  ord_af_obs <- which.max(depth_fun(dist_af_obs))
  ord_eu_obs <- which.max(depth_fun(dist_eu_obs))
  
  observed <- full_wass_matrix[
    idx_african[ord_af_obs],
    idx_europe[ord_eu_obs]
  ]
  
  # ---- Permutation distribution ----
  
  n_af <- length(idx_african)
  n_eu <- length(idx_europe)
  idx_all <- c(idx_african, idx_europe)
  
  distribution <- numeric(iter)
  
  for (b in 1:iter) {
    
    perm_idx <- sample(idx_all, replace = FALSE)
    
    africa_fake <- perm_idx[1:n_af]
    europe_fake <- perm_idx[(n_af + 1):(n_af + n_eu)]
    
    dist_af_fake <- compute_wass_submatrix(africa_fake)
    dist_eu_fake <- compute_wass_submatrix(europe_fake)
    
    ord_af <- which.max(depth_fun(dist_af_fake))
    ord_eu <- which.max(depth_fun(dist_eu_fake))
    
    distribution[b] <- full_wass_matrix[
      africa_fake[ord_af],
      europe_fake[ord_eu]
    ]
  }
  
  mean(distribution >= observed)
}



pvalue_agehist_rev_fast <- function(iter, depth_fun, n_reverse) {
  
  idx_af_rev <- idx_african
  idx_eu_rev <- idx_europe
  
  swap_idx <- sample(seq_along(idx_eu_rev), n_reverse)
  
  temp <- idx_eu_rev[swap_idx]
  idx_eu_rev[swap_idx] <- idx_af_rev[swap_idx]
  idx_af_rev[swap_idx] <- temp
  
  # ---- Observed ----
  
  dist_af_obs <- compute_wass_submatrix(idx_af_rev)
  dist_eu_obs <- compute_wass_submatrix(idx_eu_rev)
  
  ord_af_obs <- which.max(depth_fun(dist_af_obs))
  ord_eu_obs <- which.max(depth_fun(dist_eu_obs))
  
  observed <- full_wass_matrix[
    idx_af_rev[ord_af_obs],
    idx_eu_rev[ord_eu_obs]
  ]
  
  # ---- Permutation distribution ----
  
  n_af <- length(idx_af_rev)
  n_eu <- length(idx_eu_rev)
  idx_all <- c(idx_af_rev, idx_eu_rev)
  
  distribution <- numeric(iter)
  
  for (b in 1:iter) {
    
    perm_idx <- sample(idx_all, replace = FALSE)
    
    africa_fake <- perm_idx[1:n_af]
    europe_fake <- perm_idx[(n_af + 1):(n_af + n_eu)]
    
    dist_af_fake <- compute_wass_submatrix(africa_fake)
    dist_eu_fake <- compute_wass_submatrix(europe_fake)
    
    ord_af <- which.max(depth_fun(dist_af_fake))
    ord_eu <- which.max(depth_fun(dist_eu_fake))
    
    distribution[b] <- full_wass_matrix[
      africa_fake[ord_af],
      europe_fake[ord_eu]
    ]
  }
  
  mean(distribution >= observed)
}


iter <- 500
cores <- detectCores() - 1






results_df <- do.call(rbind,
                      mclapply(1:length(depth_models), function(i) {
                        
                        set.seed(2223 + i)
                        
                        data.frame(
                          Model = depth_names[i],
                          PValue = pvalue_age_hist_fast(iter, depth_models[[i]])
                        )
                        
                      }, mc.cores = cores)
)

print(results_df)


number_of_patches <- 10
iter <- 500

results_list <- mclapply(1:number_of_patches, function(batch) {
  
  set.seed(4000 + batch)
  
  res12 <- do.call(rbind,
                   lapply(1:length(depth_models), function(i) {
                     data.frame(
                       Model = depth_names[i],
                       PValue = pvalue_agehist_rev_fast(iter, depth_models[[i]], 12)
                     )
                   })
  )
  
  res16 <- do.call(rbind,
                   lapply(1:length(depth_models), function(i) {
                     data.frame(
                       Model = depth_names[i],
                       PValue = pvalue_agehist_rev_fast(iter, depth_models[[i]], 16)
                     )
                   })
  )
  
  list(res12 = res12, res16 = res16)
  
}, mc.cores = cores)



results_df_12 <- bind_rows(lapply(results_list, `[[`, "res12")) %>%
  group_by(Model) %>%
  summarise(PValue = mean(PValue))

results_df_16 <- bind_rows(lapply(results_list, `[[`, "res16")) %>%
  group_by(Model) %>%
  summarise(PValue = mean(PValue))
