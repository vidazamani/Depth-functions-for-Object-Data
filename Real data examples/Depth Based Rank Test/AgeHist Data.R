library(HistDAWass)
library(tidyr)
library(Rcpp) ## you need to have the latest version of Rcpp
library(MetricDepth)
library(ggpattern) 
library(patchwork)
library(ggpubr)
library(dplyr)
library(ggplot2)



remotes::install_github('vidazamani/Depth-functions-for-Object-Data/RelativeDepth')
library(RelativeDepth)


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

####  Wasserstein Distance Function

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


depth_names <- c("MOD3", "MOD2", "MHD", "MLD", "MSD")
depth_models <- list(MOD3_cpp, MOD2_cpp, MHD_cpp, MLD_cpp, MSD_cpp)


depth_relative_models <- list(MOD3_relative,
                              MOD2_relative,
                              MHD_relative,
                              MLD_relative,
                              MSD_relative)



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

for (i in 1:length(depth_models)) {
  pvalw <- rank_test_agehist(depth_models[[i]])$pvaluew
  pvalk <- rank_test_agehist(depth_models[[i]])$pvaluek
  results_rank <- rbind(results_rank,
                        data.frame(Model = depth_names[i],
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

for (i in 1:length(depth_models)) {
  out <- compute_K_type(depth_models[[i]])
  
  results_K <- rbind(results_K,
                     data.frame(Model = depth_names[i],
                                K_stat = out$K_statistic,
                                PValue = out$p_value))
}

print(results_K)


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


for (i in 1:length(depth_relative_models)) {
  
  out <- H_statistic_relative(
    depth_relative_fun = depth_relative_models[[i]],
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






################### Swaps #####################

compute_K_with_groups <- function(depth_fun, idx_group1, idx_group2) {
  
  idx_all <- c(idx_group1, idx_group2)
  
  n1 <- length(idx_group1)
  n2 <- length(idx_group2)
  n  <- n1 + n2
  
  dist_mat_all <- compute_wass_matrix(Age_Pyramids_2014, idx_all)
  
  depth_values <- depth_fun(dist_mat_all)
  
  ranks <- rank(depth_values, ties.method = "average")
  
  R1_sum <- sum(ranks[1:n1])
  R2_sum <- sum(ranks[(n1+1):n])
  
  K <- (12 / (n * (n + 1))) *
    (R1_sum^2 / n1 + R2_sum^2 / n2) -
    3 * (n + 1)
  
  p_value <- 1 - pchisq(K, df = 1)
  
  return(p_value)
}



compute_H_with_groups <- function(depth_relative_fun,
                                  idx_group1,
                                  idx_group2) {
  
  idx_all <- c(idx_group1, idx_group2)
  
  n1 <- length(idx_group1)
  n2 <- length(idx_group2)
  n  <- n1 + n2
  
  dist_mat_all <- compute_wass_matrix(Age_Pyramids_2014, idx_all)
  
  compute_Hk <- function(ref_positions) {
    
    depth_vals <- depth_relative_fun(dist_mat_all, ref_positions)
    
    ranks <- rank(depth_vals, ties.method = "average")
    
    R1_bar <- mean(ranks[1:n1])
    R2_bar <- mean(ranks[(n1+1):n])
    
    Hk <- (12 / (n * (n + 1))) *
      (n1 * (R1_bar - (n + 1)/2)^2 +
         n2 * (R2_bar - (n + 1)/2)^2)
    
    return(Hk)
  }
  
  H1 <- compute_Hk(1:n1)
  H2 <- compute_Hk((n1+1):n)
  
  H_value <- (H1 + H2) / 2
  
  p_value <- 1 - pchisq(H_value, df = 1)
  
  return(p_value)
}

contamination_rank_tests <- function(n_reverse,
                                     n_batches = 10,
                                     depth_fun,
                                     depth_relative_fun) {
  
  pvals_K <- c()
  pvals_H <- c()
  
  for (batch in 1:n_batches) {
    
    set.seed(4000 + batch)
    
    idx_africa_rev <- idx_african
    idx_europe_rev <- idx_europe
    
    idx_swap <- sample(length(idx_europe_rev), n_reverse)
    
    # Swap labels
    temp <- idx_europe_rev[idx_swap]
    idx_europe_rev[idx_swap] <- idx_africa_rev[idx_swap]
    idx_africa_rev[idx_swap] <- temp
    
    # K-type
    pvals_K[batch] <- compute_K_with_groups(depth_fun,
                                            idx_africa_rev,
                                            idx_europe_rev)
    
    # H-type
    pvals_H[batch] <- compute_H_with_groups(depth_relative_fun,
                                            idx_africa_rev,
                                            idx_europe_rev)
  }
  
  return(list(K_mean = mean(pvals_K),
              H_mean = mean(pvals_H)))
}






results_12 <- data.frame()

for (i in 1:length(depth_models)) {
  
  res <- contamination_rank_tests(
    n_reverse = 12,
    depth_fun = depth_models[[i]],
    depth_relative_fun = depth_relative_models[[i]]
  )
  
  results_12 <- rbind(results_12,
                      data.frame(Model = depth_names[i],
                                 K_PValue = res$K_mean,
                                 H_PValue = res$H_mean))
}

results_16 <- data.frame()

for (i in 1:length(depth_models)) {
  
  res <- contamination_rank_tests(
    n_reverse = 16,
    depth_fun = depth_models[[i]],
    depth_relative_fun = depth_relative_models[[i]]
  )
  
  results_16 <- rbind(results_16,
                      data.frame(Model = depth_names[i],
                                 K_PValue = res$K_mean,
                                 H_PValue = res$H_mean))
}


print(results_12)
print(results_16)



####### Visualization 

df_12 <- results_12 %>%
  tidyr::pivot_longer(cols = c(K_PValue, H_PValue),
                      names_to = "Test",
                      values_to = "PValue") %>%
  dplyr::mutate(
    Test = dplyr::recode(Test,
                         K_PValue = "K-type",
                         H_PValue = "H-type"),
    Scenario = "12 swaps"
  )


df_16 <- results_16 %>%
  tidyr::pivot_longer(cols = c(K_PValue, H_PValue),
                      names_to = "Test",
                      values_to = "PValue") %>%
  dplyr::mutate(
    Test = dplyr::recode(Test,
                         K_PValue = "K-type",
                         H_PValue = "H-type"),
    Scenario = "16 swaps"
  )


df_clean <- df_plot %>%
  dplyr::select(Model, Test, PValue) %>%
  dplyr::mutate(Scenario = "No swaps")


df_all <- dplyr::bind_rows(df_clean, df_12, df_16)

df_all <- df_all %>%
  dplyr::mutate(logP = -log10(PValue))

df_all$Scenario <- factor(df_all$Scenario,
                          levels = c("No swaps",
                                     "12 swaps",
                                     "16 swaps"))

p_all <- ggplot(df_all,
                aes(x = Model,
                    y = logP,
                    fill = Test)) +
  
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           color = "black",
           linewidth = 0.4) +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             linewidth = 0.8,
             color = "blue") +
  
  scale_fill_manual(values = c("K-type" = "grey40",
                               "H-type" = "grey75")) +
  
  facet_wrap(~ Scenario, nrow = 1, strip.position = "top") +
  
  labs(
    x = "",
    y = expression(-log[10](italic(p_value))),
    fill = ""
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid = element_blank()
  )


print(p_all)
