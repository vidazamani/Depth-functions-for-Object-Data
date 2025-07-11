library(HistDAWass)
library(tidyr)
library(Rcpp) ## you need to have the latest version of Rcpp
library(MetricDepth)
library(ggpattern) 
library(patchwork)
library(ggpubr)


Age_Pyramids_2014


# Pull out the M matrix of distributionH objects
M <- attr(Age_Pyramids_2014, "M")
countries <- rownames(M)

# Create a data frame with every country and its index
df_regions <- data.frame(
  country = countries,
  index   = seq_along(countries),
  stringsAsFactors = FALSE
)


# 2. Define Factbook subregion vectors
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

###### Examples shwoing how to use the function above 

# Distance matrix for each group
# dist_mat_eastern <- compute_wass_matrix(Age_Pyramids_2014,idx_east)
# dist_mat_western <- compute_wass_matrix(Age_Pyramids_2014,idx_west)
# 
# 
# # Deepest object find by depths functions
# ordes <- order(MOD3_cpp(dist_mat_eastern),decreasing = TRUE)[1]
# deepestES <- df_regions[df_regions$index == idx_east[ordes],1]
# deepestEShist <- Age_Pyramids_2014[idx_east[ordes]]@M[[1]]
# 
# ordws <- order(MSD_cpp(dist_mat_western),decreasing = TRUE)[1]
# deepestWS <- df_regions[df_regions$index == idx_west[ordws],1]
# deepestWShist <- Age_Pyramids_2014[idx_west[ordws]]@M[[1]]





pvalue_age_hist <- function(iter, depth_fun) {
  
  distribution=c()
  
  for (i in 1:iter) {
    
    
    ###### Western vs Eastern
    # 
    # fake_labels <- sample(c(rep("es",16),rep("ws",24)), 40, replace = FALSE)
    # 
    # east_fake <- df_europe$index[fake_labels %in% c("es")]
    # 
    # west_fake <- df_europe$index[fake_labels %in% c("ws")]
    
    
    ###### Africa vs Europe
    
    fake_labels <- sample(c(rep("af",55),rep("eu",40)), 95, replace = FALSE)
    
    df_all <- rbind(df_african,df_europe)
    
    africa_fake <- df_all$index[fake_labels %in% c("af")]
    
    europe_fake <- df_all$index[fake_labels %in% c("eu")]
    
    
    ###### Computing Distance for fake indices 
    
    ##### Western vs Eastern
    
    # dist_mat_eastern_fake <- compute_wass_matrix(Age_Pyramids_2014,na.omit(east_fake))
    # dist_mat_western_fake <- compute_wass_matrix(Age_Pyramids_2014,na.omit(west_fake))
    
    ###### Africa vs Europe
    
    dist_mat_africa_fake <- compute_wass_matrix(Age_Pyramids_2014,na.omit(africa_fake))
    dist_mat_europe_fake <- compute_wass_matrix(Age_Pyramids_2014,na.omit(europe_fake))
    
    
    ##################### Difference in means
    
    #### Western vs Eastern
    
    # ordes <- order(depth_fun(dist_mat_eastern_fake),decreasing = TRUE)[1]
    # deepestEShist <- Age_Pyramids_2014[idx_east[ordes]]@M[[1]]
    # 
    # 
    # 
    # ordws <- order(depth_fun(dist_mat_western_fake),decreasing = TRUE)[1]
    # deepestWShist <- Age_Pyramids_2014[idx_west[ordws]]@M[[1]]
    # 
    # 
    # 
    # distribution[i] <- WassSqDistH(deepestEShist,deepestWShist)
    # 
    
    ###### Africa vs Europe
    
    ordaf <- order(depth_fun(dist_mat_africa_fake),decreasing = TRUE)[1]
    deepestAFhist <- Age_Pyramids_2014[africa_fake[ordaf]]@M[[1]]
    
    
    
    ordeu <- order(depth_fun(dist_mat_europe_fake),decreasing = TRUE)[1]
    deepestEUhist <- Age_Pyramids_2014[europe_fake[ordeu]]@M[[1]]
    
    
    
    distribution[i] <- sqrt(WassSqDistH(deepestAFhist,deepestEUhist))
    
    
  }
  
  #Difference in means for real observed data
  
  
  #### Western vs Eastern
  
  # dist_mat_eastern <- compute_wass_matrix(Age_Pyramids_2014,idx_east)
  # dist_mat_western <- compute_wass_matrix(Age_Pyramids_2014,idx_west)
  # 
  
  
  # ordes_r <- order(depth_fun(dist_mat_eastern),decreasing = TRUE)[1]
  # 
  # ordws_r <- order(depth_fun(dist_mat_western),decreasing = TRUE)[1]
  # 
  # observed <- WassSqDistH(Age_Pyramids_2014[idx_east[ordes_r]]@M[[1]],
  #                         Age_Pyramids_2014[idx_west[ordws_r]]@M[[1]])
  # 
  # 
  
  
  #### Africa vs Europe 
  
  dist_mat_africa <- compute_wass_matrix(Age_Pyramids_2014,idx_african)
  dist_mat_europe <- compute_wass_matrix(Age_Pyramids_2014,idx_europe)
  
  
  ordaf_r <- order(depth_fun(dist_mat_africa),decreasing = TRUE)[1]
  
  ordeu_r <- order(depth_fun(dist_mat_europe),decreasing = TRUE)[1]
  
  observed <- sqrt(WassSqDistH(Age_Pyramids_2014[idx_african[ordaf_r]]@M[[1]],
                               Age_Pyramids_2014[idx_europe[ordeu_r]]@M[[1]]))
  
  
  
  pvalue <- sum(abs(distribution) >= abs(c(observed)))/(iter)
  
  return(pvalue)
  
}


start_time <- Sys.time()
set.seed(2223)


##### Parameters 

iter = 10

model_names <- c("MOD3", "MOD2", "MHD", "MLD", "MSD")
models <- list(MOD3_cpp, MOD2_cpp, MHD_cpp, MLD_cpp, MSD_cpp)


combinations <- expand.grid(ModelName = model_names,
                            stringsAsFactors = FALSE)


# Apply the function for each combination and store the results in a data frame
results_df <- do.call(rbind, lapply(1:nrow(combinations), function(i) {
  model <- models[[which(model_names == combinations$ModelName[i])]]
  data.frame(Model = combinations$ModelName[i], PValue = pvalue_age_hist(iter, model))
}))

# Print the final results
print(results_df)


end_time <- Sys.time()

running_time <- end_time-start_time











#################### Reverse #######################


pvalue_agehist_rev <- function(iter, depth_fun, n_reverse) {
  
  #### Western vs Eastern
  
  # df_eastern_rev <- df_eastern$country
  # df_western_rev <- df_western$country
  # 
  # 
  # idx_east_rev <- sample(length(df_eastern_rev), n_reverse)
  # 
  # # Perform the swap
  # temp <- df_eastern_rev[idx_east_rev]
  # df_eastern_rev[idx_east_rev] <- df_western_rev[idx_east_rev]
  # df_western_rev[idx_east_rev] <- temp
  # 
  # 
  # df_es_rev <- df_eastern
  # df_es_rev$country = df_eastern_rev
  # 
  # 
  # df_ws_rev <- df_western
  # df_ws_rev$country <- df_western_rev
  # 
  # #### Africa vs Europe
  # 
  df_africa_rev <- df_african$country
  df_europe_rev <- df_europe$country
  
  idx_europe_rev <- idx_europe
  idx_african_rev <- idx_african
  
  idx_swap <- sample(length(df_europe_rev), n_reverse)
  
  # Perform the swap
  temp <- idx_europe[idx_swap]
  idx_europe_rev[idx_swap] <- idx_african_rev[idx_swap]
  idx_african_rev[idx_swap] <- temp
  
  
  df_all_now <- rbind(df_african, df_europe)
  
  df_af_rev <- data.frame()
  for(i in 1:length(idx_african_rev)){
    df_af_rev <- bind_rows(df_af_rev, df_all_now[df_all_now$index == idx_african_rev[i], ])
  }
  df_af_rev
  
  df_eu_rev <- data.frame()
  for(i in 1:length(idx_europe_rev)){
    df_eu_rev <- bind_rows(df_eu_rev, df_all_now[df_all_now$index == idx_europe_rev[i], ])
  }
  df_eu_rev
  
  
  distribution <- c()
  
  
  for (i in 1:iter) {
    
    #### Western vs Eastern 
    
    # fake_labels <- sample(c(rep("es",length(df_eastern_rev)),
    #                         rep("ws",length(df_western_rev))), 40, replace = FALSE)
    # 
    # 
    # east_fake <- df_europe$index[fake_labels %in% c("es")]
    # 
    # west_fake <- df_europe$index[fake_labels %in% c("ws")]
    # 
    
    #### Africa vs Europe
    
    fake_labels <- sample(c(rep("af",nrow(df_af_rev)),
                            rep("eu",nrow(df_eu_rev))), 95, replace = FALSE)
    
    df_all <- rbind(df_af_rev, df_eu_rev)
    
    africa_fake <- df_all$index[fake_labels %in% c("af")]
    
    europe_fake <- df_all$index[fake_labels %in% c("eu")]
    
    
    ###### Computing Distance for fake indices 
    
    #### Western vs Eeastern
    
    # dist_mat_eastern_fake <- compute_wass_matrix(Age_Pyramids_2014,na.omit(east_fake))
    # dist_mat_western_fake <- compute_wass_matrix(Age_Pyramids_2014,na.omit(west_fake))
    # 
    
    #### Africa vs Europe
    
    dist_mat_africa_fake <- compute_wass_matrix(Age_Pyramids_2014,na.omit(africa_fake))
    dist_mat_europe_fake <- compute_wass_matrix(Age_Pyramids_2014,na.omit(europe_fake))
    
    
    ################### Difference in means
    
    
    #### Western vs Eastern 
    
    # ordes <- order(depth_fun(dist_mat_eastern_fake),decreasing = TRUE)[1]
    # deepestEShist <- Age_Pyramids_2014[idx_east[ordes]]@M[[1]]
    # 
    # 
    # 
    # ordws <- order(depth_fun(dist_mat_western_fake),decreasing = TRUE)[1]
    # deepestWShist <- Age_Pyramids_2014[idx_west[ordws]]@M[[1]]
    # 
    # 
    # 
    # distribution[i] <- WassSqDistH(deepestEShist,deepestWShist)
    # 
    
    
    #### Africa vs Europe 
    
    ordaf <- order(depth_fun(dist_mat_africa_fake),decreasing = TRUE)[1]
    deepestAFhist <- Age_Pyramids_2014[africa_fake[ordaf]]@M[[1]]
    
    
    
    ordeu <- order(depth_fun(dist_mat_europe_fake),decreasing = TRUE)[1]
    deepestEUhist <- Age_Pyramids_2014[europe_fake[ordeu]]@M[[1]]
    
    
    
    distribution[i] <- sqrt(WassSqDistH(deepestAFhist,deepestEUhist))
    
    
    
  }
  
  
  
  ##### Difference in means for real observed data
  
  #### Western vs Eastern
  
  # dist_mat_eastern <- compute_wass_matrix(Age_Pyramids_2014,df_es_rev$index)
  # dist_mat_western <- compute_wass_matrix(Age_Pyramids_2014,df_ws_rev$index)
  # 
  # 
  # 
  # ordes_r <- order(depth_fun(dist_mat_eastern),decreasing = TRUE)[1]
  # 
  # ordws_r <- order(depth_fun(dist_mat_western),decreasing = TRUE)[1]
  # 
  # observed <- WassSqDistH(Age_Pyramids_2014[df_es_rev$index[ordes_r]]@M[[1]],
  #                         Age_Pyramids_2014[df_ws_rev$index[ordws_r]]@M[[1]])
  # 
  
  #### Africa vs Europe
  
  dist_mat_africa <- compute_wass_matrix(Age_Pyramids_2014,df_af_rev$index)
  dist_mat_europe <- compute_wass_matrix(Age_Pyramids_2014,df_eu_rev$index)
  
  
  
  ordaf_r <- order(depth_fun(dist_mat_africa),decreasing = TRUE)[1]
  
  ordeu_r <- order(depth_fun(dist_mat_europe),decreasing = TRUE)[1]
  
  observed <- sqrt(WassSqDistH(Age_Pyramids_2014[df_af_rev$index[ordaf_r]]@M[[1]],
                               Age_Pyramids_2014[df_eu_rev$index[ordeu_r]]@M[[1]]))
  
  
  pvalue <- sum(abs(distribution) >= abs(c(observed)))/(iter)
  
  return(pvalue)
  
}



number_of_patches <- 10
iter <- 100

results_df_7_final <- data.frame()
results_df_12_final <- data.frame()


for(batch in 1:number_of_patches){
  start_time <- Sys.time()
  set.seed(2222 + batch)
  
  
  ##### Parameters for 7 swaps
  
  
  n_reverse <- 12
  
  
  model_names <- c("MOD3", "MOD2", "MHD", "MLD", "MSD")
  models <- list(MOD3_cpp, MOD2_cpp, MHD_cpp, MLD_cpp, MSD_cpp)
  
  
  combinations <- expand.grid(ModelName = model_names,
                              stringsAsFactors = FALSE)
  
  
  # Apply the function for each combination and store the results in a data frame
  results_df_7 <- do.call(rbind, lapply(1:nrow(combinations), function(i) {
    model <- models[[which(model_names == combinations$ModelName[i])]]
    data.frame(Model = combinations$ModelName[i],
               PValue = pvalue_agehist_rev(iter, model,n_reverse))
  }))
  
  # Print the final results
  # print(results_df_7)
  
  results_df_7_final <- bind_rows(results_df_7_final, results_df_7)
  
  end_time <- Sys.time()
  
  running_time <- end_time-start_time
  
  
  ##### Parameters for 15 swaps
  
  n_reverse <- 16
  
  
  model_names <- c("MOD3", "MOD2", "MHD", "MLD", "MSD")
  models <- list(MOD3_cpp, MOD2_cpp, MHD_cpp, MLD_cpp, MSD_cpp)
  
  
  combinations <- expand.grid(ModelName = model_names,
                              stringsAsFactors = FALSE)
  
  
  # Apply the function for each combination and store the results in a data frame
  results_df_12 <- do.call(rbind, lapply(1:nrow(combinations), function(i) {
    model <- models[[which(model_names == combinations$ModelName[i])]]
    data.frame(Model = combinations$ModelName[i],
               PValue = pvalue_agehist_rev(iter, model,n_reverse))
  }))
  
  # Print the final results
  # print(results_df_12)
  
  results_df_12_final <- bind_rows(results_df_12_final, results_df_12)
  
  end_time <- Sys.time()
  
  running_time <- end_time-start_time
  
  
  print(batch)
}



results_df_7 <- results_df_7_final %>%
  group_by(Model) %>%
  summarise(PValue = mean(PValue))

results_df_12 <- results_df_12_final %>%
  group_by(Model) %>%
  summarise(PValue = mean(PValue))



#



#


My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 8),
  axis.title.y = element_text(size = 12))



p1 = ggplot(results_df, aes(x = Model, y = PValue, fill = Model)) +
  geom_bar_pattern(stat = "identity", position = "dodge", color = "black",
                   pattern_density = 0.4,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02
  ) +
  scale_fill_grey(start = 0.1, end = 1) +  
  scale_pattern_manual(values = c("none", "crosshatch", "stripe", "circle", "none")) +  # Custom patterns
  labs(title = "",
       x = expression(paste('Without swaps')),
       y = expression(paste(italic(p),'-',values)), fill = " ", pattern = " ") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),
        legend.text=element_text(size=11),
        strip.text.x = element_text(size = 25, 
                                    colour = "black",
                                    angle = 90),
        strip.text.y = element_text(size = 25, 
                                    colour = "black",
                                    angle = 90),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE),
         pattern = guide_legend(nrow = 1, byrow = TRUE))+
  geom_hline(yintercept=0.05, linetype="dashed", color = "blue",
             linewidth=1)+ My_Theme+ 
  scale_x_discrete(labels = c('0' = expression(infinity)))


p2 = ggplot(results_df_7, aes(x = Model, y = PValue, fill = Model)) +
  geom_bar_pattern(stat = "identity", position = "dodge", color = "black",
                   pattern_density = 0.4,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02
  ) +
  scale_fill_grey(start = 0.1, end = 1) +  
  scale_pattern_manual(values = c("none", "crosshatch", "stripe", "circle", "none")) +  # Custom patterns
  labs(title = "",
       x = expression(paste('12 swaps')),
       y = expression(paste(italic(p),'-',values)), fill = " ", pattern = " ") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),
        legend.text=element_text(size=11),
        strip.text.x = element_text(size = 25, 
                                    colour = "black",
                                    angle = 90),
        strip.text.y = element_text(size = 25, 
                                    colour = "black",
                                    angle = 90),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE),
         pattern = guide_legend(nrow = 1, byrow = TRUE))+
  geom_hline(yintercept=0.05, linetype="dashed", color = "blue",
             linewidth=1)+ My_Theme+ 
  scale_x_discrete(labels = c('0' = expression(infinity)))



p3 = ggplot(results_df_12, aes(x = Model, y = PValue, fill = Model)) +
  geom_bar_pattern(stat = "identity", position = "dodge", color = "black",
                   pattern_density = 0.4,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02
  ) +
  scale_fill_grey(start = 0.1, end = 1) +  
  scale_pattern_manual(values = c("none", "crosshatch", "stripe", "circle", "none")) +  # Custom patterns
  labs(title = "",
       x = expression(paste('16 swaps')),
       y = expression(paste(italic(p),'-',values)), fill = " ", pattern = " ") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),
        legend.text=element_text(size=11),
        strip.text.x = element_text(size = 25, 
                                    colour = "black",
                                    angle = 90),
        strip.text.y = element_text(size = 25, 
                                    colour = "black",
                                    angle = 90),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE),
         pattern = guide_legend(nrow = 1, byrow = TRUE))+
  geom_hline(yintercept=0.05, linetype="dashed", color = "blue",
             linewidth=1)+ My_Theme+ 
  scale_x_discrete(labels = c('0' = expression(infinity)))




pdf("/Users/jomivi/Library/Mobile Documents/com~apple~CloudDocs/Work/Supervision/Vida/Paper_1_revision/reversed and original.pdf", 10 ,6)
ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="none")
dev.off()
