library(HistDAWass)
library(tidyr)
library(Rcpp) ## you need to have the latest version of Rcpp
library(MetricDepth)
library(ggpattern) 
library(patchwork)
library(ggpubr)
library(dplyr)
library(ggplot2)





# Pull out the M matrix of distributionH objects
M <- attr(Age_Pyramids_2014, "M")
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

#### Put Eastern_index or Western_index as a group input

compute_wass_matrix <- function(hist_data,group) {
  n <- length(group) 
  dist_matrix <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      d <- WassSqDistH(hist_data[group[i]]@M[[1]],
                       hist_data[group[j]]@M[[1]])
      dist_matrix[i, j] <- d
      dist_matrix[j, i] <- d
    }
  }
  return(sqrt(dist_matrix))
}



##### Parameters 


model_names <- c("MOD3", "MOD2", "MHD", "MLD", "MSD")
models <- list(MOD3_cpp, MOD2_cpp, MHD_cpp, MLD_cpp, MSD_cpp)





################### Rank Test #################

rank_test_agehist <- function(depth_fun) {
  
  idx_all <- c(idx_african, idx_europe)
  
  dist_mat_all <- compute_wass_matrix(Age_Pyramids_2014, idx_all)
  
  depth_values <- depth_fun(dist_mat_all)
  
  group_labels <- c(rep("Africa", length(idx_african)),
                    rep("Europe", length(idx_europe)))
  
  
  test <- wilcox.test(depth_values ~ group_labels)
  test2 <- kruskal.test(depth_values ~ group_labels)
  
  ## or
  # test <- wilcox.test(depth_values[group_labels == "Africa"],
  #                     depth_values[group_labels == "Europe"],
  #                     alternative = "two.sided")
  
  output <- list("pvaluew" = test$p.value, 
                 "pvaluek" = test2$p.value,
                 "depths" = depth_values, "grouplab" = group_labels)
  return(output)
}




results_rank <- data.frame()

for (i in 1:length(models)) {
  pvalw <- rank_test_agehist(models[[i]])$pvaluew
  pvalk <- rank_test_agehist(models[[i]])$pvaluek
  results_rank <- rbind(results_rank,
                        data.frame(Model = model_names[i],
                                   PValuew = pvalw, PValuek = pvalk))
}

print(results_rank)

########### K type statistics in Chenouri paper ##############################################
## it is completely similar to kruskal.test when we have 2 samples
## Under large samples: K and H become asymptotically equivalent.
## kruskal test becomes wilcoxon test when we have 2 samples


compute_K_type <- function(depth_fun) {
  
  # Combine indices
  idx_all <- c(idx_african, idx_europe)
  
  n1 <- length(idx_african)
  n2 <- length(idx_europe)
  n  <- n1 + n2
  
  # Step 1: compute pooled distance matrix
  dist_mat_all <- compute_wass_matrix(Age_Pyramids_2014, idx_all)
  
  # Step 2: compute pooled depths
  depth_values <- depth_fun(dist_mat_all)
  
  # Step 3: rank depths
  ranks <- rank(depth_values)
  
  # Step 4: rank sums
  R1_sum <- sum(ranks[1:n1])
  R2_sum <- sum(ranks[(n1+1):n])
  
  # Step 5: compute K statistic
  K <- (12 / (n * (n + 1))) * 
    (R1_sum^2 / n1 + R2_sum^2 / n2) -
    3 * (n + 1)
  
  # Asymptotic p-value
  p_value <- 1 - pchisq(K, df = 1)
  
  return(list(K_statistic = K, p_value = p_value))
}

results_K <- data.frame()

for (i in 1:length(models)) {
  out <- compute_K_type(models[[i]])
  
  results_K <- rbind(results_K,
                     data.frame(Model = model_names[i],
                                K_stat = out$K_statistic,
                                PValue = out$p_value))
}

print(results_K)

################# Relative Depth Functions ############################

MLD_relative <- function(D, ref_idx) {
  
  n_total <- nrow(D)
  n_ref   <- length(ref_idx)
  
  lens_depth <- rep(0, n_total)
  
  for (p in 1:n_total) {
    
    s <- 0
    
    for (a in 1:(n_ref-1)) {
      for (b in (a+1):n_ref) {
        
        i <- ref_idx[a]
        j <- ref_idx[b]
        
        maximum <- max(D[i,p], D[j,p])
        
        if (D[i,j] > maximum + 1e-6) {
          s <- s + 1
        }
      }
    }
    
    lens_depth[p] <- (1/choose(n_ref,2)) * s
  }
  
  return(lens_depth)
}


MLD_relative <- function(D, ref_idx) {
  
  n_total <- nrow(D)
  n_ref   <- length(ref_idx)
  
  lens_depth <- rep(0, n_total)
  
  for (p in 1:n_total) {
    
    s <- 0
    
    for (a in 1:(n_ref-1)) {
      for (b in (a+1):n_ref) {
        
        i <- ref_idx[a]
        j <- ref_idx[b]
        
        maximum <- max(D[i,p], D[j,p])
        
        if (D[i,j] > maximum + 1e-6) {
          s <- s + 1
        }
      }
    }
    
    lens_depth[p] <- (1/choose(n_ref,2)) * s
  }
  
  return(lens_depth)
}


MHD_relative <- function(D, ref_idx) {
  
  n_total <- nrow(D)
  n_ref   <- length(ref_idx)
  
  p_mat <- matrix(0, n_ref, n_ref)
  
  # Build halfspaces using reference group only
  for (a in 1:(n_ref-1)) {
    for (b in (a+1):n_ref) {
      
      i <- ref_idx[a]
      j <- ref_idx[b]
      
      s <- 0
      
      for (k in ref_idx) {
        if (D[k,i] <= D[k,j] + 1e-6) {
          s <- s + 1
        }
      }
      
      p_mat[a,b] <- s / n_ref
    }
  }
  
  depth_vals <- rep(0, n_total)
  
  for (y in 1:n_total) {
    
    Q <- c()
    
    for (a in 1:(n_ref-1)) {
      for (b in (a+1):n_ref) {
        
        i <- ref_idx[a]
        j <- ref_idx[b]
        
        if (D[y,i] <= D[y,j] + 1e-6) {
          Q <- c(Q, p_mat[a,b])
        }
      }
    }
    
    depth_vals[y] <- min(Q)
  }
  
  return(depth_vals)
}


MSD_relative <- function(D, ref_idx) {
  
  n_total <- nrow(D)
  n_ref   <- length(ref_idx)
  
  res <- rep(0, n_total)
  
  for (k in 1:n_total) {
    
    res_now <- 0
    
    for (i in ref_idx) {
      for (j in ref_idx) {
        
        if (i != j && D[k,i] > 1e-6 && D[k,j] > 1e-6) {
          
          temp <- D[k,i] / D[k,j]
          
          res_now <- res_now + 
            temp + 1/temp - D[i,j]^2 / (D[k,i] * D[k,j])
        }
      }
    }
    
    res[k] <- res_now
  }
  
  res <- res / (n_ref^2)
  
  return(1 - 0.5 * res)
}


MOD2_relative <- function(D, ref_idx) {
  
  n_total <- nrow(D)
  n_ref   <- length(ref_idx)
  
  ojadepth <- rep(0, n_total)
  
  for (w in 1:n_total) {
    
    area <- 0
    
    for (a in 1:(n_ref-1)) {
      for (b in (a+1):n_ref) {
        
        i <- ref_idx[a]
        j <- ref_idx[b]
        
        S <- matrix(c(
          (D[i,w])^2,
          -0.5*((D[j,i])^2 - D[j,w]^2 - D[i,w]^2),
          -0.5*((D[i,j])^2 - D[i,w]^2 - D[j,w]^2),
          (D[j,w])^2
        ), 2, 2)
        
        detS <- det(S)
        if (detS > 1e-6)
          area <- area + sqrt(detS)
      }
    }
    
    ojadepth[w] <- 1 / (1 + (1/(0.5*(n_ref^2 - n_ref))) * area)
  }
  
  return(ojadepth)
}

MOD3_relative <- function(D, ref_idx) {
  
  n_total <- nrow(D)
  n_ref   <- length(ref_idx)
  
  ojadepth <- rep(0, n_total)
  
  for (w in 1:n_total) {
    
    area <- 0
    
    for (a in 1:(n_ref-2)) {
      for (b in (a+1):(n_ref-1)) {
        for (c in (b+1):n_ref) {
          
          i <- ref_idx[a]
          j <- ref_idx[b]
          k <- ref_idx[c]
          
          S <- matrix(c(
            (D[i,w])^2,
            -0.5*((D[j,i])^2 - D[j,w]^2 - D[i,w]^2),
            -0.5*((D[k,i])^2 - D[k,w]^2 - D[i,w]^2),
            -0.5*((D[i,j])^2 - D[i,w]^2 - D[j,w]^2),
            (D[j,w])^2,
            -0.5*((D[k,j])^2 - D[k,w]^2 - D[j,w]^2),
            -0.5*((D[i,k])^2 - D[i,w]^2 - D[k,w]^2),
            -0.5*((D[j,k])^2 - D[j,w]^2 - D[k,w]^2),
            (D[k,w])^2
          ), 3, 3)
          
          detS <- det(S)
          
          if (detS > 1e-6)
            area <- area + sqrt(detS)
        }
      }
    }
    
    ojadepth[w] <- 1 / (1 + (1/(0.5*(n_ref-1)*(n_ref^2))) * area)
  }
  
  return(ojadepth)
}

########################################################################

# idx_all <- c(idx_african, idx_europe)
# dist_mat_all <- compute_wass_matrix(Age_Pyramids_2014, idx_all)

# n1 <- length(idx_african)
# n2 <- length(idx_europe)
# 
# ref_africa <- 1:n1
# ref_europe <- (n1+1):(n1+n2)


### OR

# ref_africa <- match(idx_african, idx_all)
# ref_europe <- match(idx_europe, idx_all)

# depth_africa <- MHD_relative(dist_mat_all, ref_idx = ref_africa)
# depth_europe <- MHD_relative(dist_mat_all, ref_idx = ref_europe)






H_statistic_relative <- function(depth_relative_fun,
                                 idx_group1,
                                 idx_group2,
                                 data_object) {
  
  # --- Step 1: Pool indices ---
  idx_all <- c(idx_group1, idx_group2)
  
  # --- Step 2: Compute full distance matrix ---
  dist_mat_all <- compute_wass_matrix(data_object, idx_all)
  
  # --- Step 3: Sample sizes ---
  n1 <- length(idx_group1)
  n2 <- length(idx_group2)
  n  <- n1 + n2
  t  <- 2
  
  # --- Step 4: Reference positions inside pooled matrix ---
  ref1 <- 1:n1
  ref2 <- (n1 + 1):(n1 + n2)
  
  # ----------------------------------------------------------
  # Function to compute H(k) for one reference group
  # ----------------------------------------------------------
  compute_Hk <- function(ref_positions) {
    
    # Compute depth relative to reference group
    depth_vals <- depth_relative_fun(dist_mat_all, ref_positions)
    
    # Rank depths
    ranks <- rank(depth_vals, ties.method = "average")
    
    # Rank sums
    R1_sum <- sum(ranks[1:n1])
    R2_sum <- sum(ranks[(n1+1):n])
    
    # Chenouri & Small formula
    Hk <- (12 / (n * (n + 1))) *
      ((R1_sum^2)/n1 + (R2_sum^2)/n2) -
      3 * (n + 1)
    
    return(Hk)
  }
  
  # --- Step 5: Compute H(1) and H(2) ---
  H1 <- compute_Hk(ref1)
  H2 <- compute_Hk(ref2)
  
  # --- Step 6: Final H statistic ---
  H_value <- (H1 + H2) / t
  
  # --- Step 7: Asymptotic p-value ---
  p_value <- 1 - pchisq(H_value, df = t - 1)
  
  return(list(H_stat = H_value,
              PValue = p_value,
              H1 = H1,
              H2 = H2))
}







results_H <- data.frame()

depth_functions <- list(
  MLD_relative,
  MHD_relative,
  MSD_relative,
  MOD2_relative,
  MOD3_relative
)

depth_names <- c("MLD","MHD","MSD","MOD2","MOD3")

for (i in 1:length(depth_functions)) {
  
  out <- H_statistic_relative(
    depth_relative_fun = depth_functions[[i]],
    idx_group1 = idx_african,
    idx_group2 = idx_europe,
    data_object = Age_Pyramids_2014
  )
  
  results_H <- rbind(results_H,
                     data.frame(Model = depth_names[i],
                                H_stat = out$H_stat,
                                PValue = out$PValue))
}

print(results_H)

################## Visualization ##############

df_K <- results_K %>%
  dplyr::select(Model, PValue) %>%
  dplyr::mutate(Test = "K-type")

df_H <- results_H %>%
  dplyr::select(Model, PValue) %>%
  dplyr::mutate(Test = "H-type")

df_plot <- dplyr::bind_rows(df_K, df_H)

df_plot <- df_plot %>%
  dplyr::mutate(logP = -log10(PValue))


p <- ggplot(df_plot, aes(x = Model, y = logP, fill = Test)) +
  
  geom_col(position = position_dodge(width = 0.65),
           width = 0.55,
           color = "black",
           linewidth = 0.4) +
  
  # Significance threshold line (alpha = 0.05)
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             linewidth = 0.6,
             color = "black") +
  
  scale_fill_grey(start = 0.35, end = 0.75) +
  
  labs(
    x = "",
    y = expression(-log[10](italic(p_value))),
    fill = "Test"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  )

print(p)




################## pseudo-H Statistics #################################
## Important note: you are still giving depth_fun a square matrix of size n × n.
## You have only Reordered rows and columns. You have NOT changed the reference distribution.
## So mathematically: D(x_i, F_k) is not being computed.
## We are still computing D(x_i, F_pooled)

# rank_test_agehist_H <- function(depth_fun) {
#   
#   idx_all <- c(idx_african, idx_europe)
#   n1 <- length(idx_african)
#   n2 <- length(idx_europe)
#   n  <- n1 + n2
#   t  <- 2
# 
#   # Precompute full distance matrix
#   dist_mat_all <- compute_wass_matrix(Age_Pyramids_2014, idx_all)
#   
#   compute_Hk <- function(ref_idx) {
#     
#     # Positions of reference group inside pooled sample
#     ref_pos <- match(ref_idx, idx_all)
#     
#     # Build a new square matrix:
#     # distances among ALL observations
#     # but depths computed relative to reference group
#     
#     # For distance-based depths, the trick is:
#     # reorder matrix so reference group is first
#     
#     perm_order <- c(ref_pos, setdiff(1:n, ref_pos))
#     dist_perm  <- dist_mat_all[perm_order, perm_order]
#     
#     # Compute depth on permuted matrix
#     depth_perm <- depth_fun(dist_perm)
#     
#     # Extract depth values corresponding to original ordering
#     depth_all <- depth_perm[order(perm_order)]
#     
#     # Rank
#     ranks <- rank(depth_all, ties.method = "average")
#     
#     # Group means
#     R1_bar <- mean(ranks[1:n1])
#     R2_bar <- mean(ranks[(n1+1):n])
#     
#     Hk <- (12 / (n * (n + 1))) * (
#       n1 * (R1_bar - (n + 1)/2)^2 +
#         n2 * (R2_bar - (n + 1)/2)^2
#     )
#     
#     return(Hk)
#   }
#   
#   H1 <- compute_Hk(idx_african)
#   H2 <- compute_Hk(idx_europe)
#   
#   H_value <- (H1 + H2) / t
#   
#   p_value <- 1 - pchisq(H_value, df = 1)
#   
#   return(list(H_stat = H_value, PValue = p_value))
# }
# 
# 
# 
# results_H <- data.frame()
# 
# for (i in 1:length(models)) {
#   out <- rank_test_agehist_H(models[[i]])
#   
#   results_H <- rbind(results_H,
#                      data.frame(Model = model_names[i],
#                                 H_stat = out$H_stat,
#                                 PValue = out$PValue))
# }
# 
# print(results_H)
# 
