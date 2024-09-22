#### Real example with Canadian weather data ####


##### All needed source and libraries ########
library(fda)
library(fda.usc) ### For computing distance
library(RColorBrewer)


source("C:/Users/vizama/Documents/1st paper/Box_Parallel/Codes/All Metric depth functions.R")


###################################################################################################
########## In this part we just want to show why this data can be divided to two groups ###########
# The full temperature data set from the Canadian weather data example
y <- CanadianWeather$dailyAv[, , 1]
my_basis <- create.fourier.basis(c(0, 365), 100)
ft <- Data2fd(1:365, y, my_basis)
ft_s <- smooth.fd(ft, fdPar(ft, 2, 1e4))

plot(ft_s,xlab = "Time (days)", ylab = "Average temperature",
     main = "Canadian Weather Data")



# Each of the weather stations belongs to one of four regions:
# Atlantic, Continental, Pacific, Arctic

region <- CanadianWeather$region



#Divide the stations into two groups ("Atlantic and Pacific" vs. "Continental and Arctic") and compute the mean and +-2 SD curves separately for the two groups.





g1rev <- region %in% c("Atlantic", "Pacific")
g2rev <- region %in% c("Continental", "Arctic")


#n_reverse: Number of elements to exchange
n_reverse = 7

indices_g1rev <- sample(length(g1rev), n_reverse)

# Perform the swap
temp <- g1rev[indices_g1rev]
g1rev[indices_g1rev] <- g2rev[indices_g1rev]
g2rev[indices_g1rev] <- temp


ft_sg1 <- ft_s[g1rev]
ft_sg2 <- ft_s[g2rev]


merged_coefs <- cbind(ft_sg1$coefs, ft_sg2$coefs)
merged_fd <- fd(coef = merged_coefs, basisobj = my_basis)






pvalue_canada_weather_rev <- function(lp_value, iter, depth_fun) {
  
  distribution <- c()

  
  for (i in 1:iter) {
    
    fake_labels <- sample(c(rep("ap",sum(g1rev)),rep("ca",sum(g2rev))), 35,
                          replace = FALSE)

    
    ap_fake <- merged_fd[fake_labels %in% c("ap")]
    
    ca_fake <- merged_fd[fake_labels %in% c("ca")]
    
    
    ###### Computing L_2 Distance 
    Dismatrixap2 <- metric.lp(fdata(ap_fake), lp = lp_value)
    
    Dismatrixca2 <- metric.lp(fdata(ca_fake), lp = lp_value)
    
    
    #Difference in means
    
    ordap <- order(depth_fun(Dismatrixap2),decreasing = TRUE)[1]
    deepestap <- fdata(ap_fake)[ordap]
    
    ordca <- order(depth_fun(Dismatrixca2),decreasing = TRUE)[1]
    deepestca <- fdata(ca_fake)[ordca]
    
    
    distribution[i] <- metric.lp(deepestap,deepestca, lp = lp_value)
    
    
  }
  
  #Difference in means for real observed data

  
  Dismatrixap_real <- metric.lp(fdata(ft_sg1), lp = lp_value)
  
  Dismatrixca_real <- metric.lp(fdata(ft_sg2), lp = lp_value)
  
  
  ordap_r <- order(depth_fun(Dismatrixap_real),decreasing = TRUE)[1]
  deepestap_real <- fdata(ft_sg1)[ordap_r]
  
  ordca_r <- order(depth_fun(Dismatrixca_real),decreasing = TRUE)[1]
  deepestca_real <- fdata(ft_sg2)[ordca_r]
  
  observed <- metric.lp(deepestap_real,deepestca_real, lp = lp_value)
  
  
  pvalue <- sum(abs(distribution) >= abs(c(observed)))/(iter)
  
  return(pvalue)
  
}

start_time <- Sys.time()
set.seed(2222)


##### Parameters 

iter = 5000
lp_values <- c(0, 1, 2, 3, 4, 5, 10, 20)
model_names <- c("MOD3", "MOD2", "MHD", "MLD", "MSD")
models <- list(MOD3, MOD2, MHD, MLD, MSD)





combinations <- expand.grid(ModelName = model_names,
                            lpvalue = lp_values,
                            stringsAsFactors = FALSE)


# Apply the function for each combination and store the results in a data frame
results_df <- do.call(rbind, lapply(1:nrow(combinations), function(i) {
  model <- models[[which(model_names == combinations$ModelName[i])]]
  lp_val <- combinations$lpvalue[i]
  data.frame(Model = combinations$ModelName[i], lpvalue = lp_val, PValue = pvalue_canada_weather_rev(lp_val, iter, model))
}))

# Print the final results
print(results_df)


end_time <- Sys.time()

running_time <- end_time-start_time

library(tidyr)

# Create the pivot table
pivot_table <- results_df %>%
  pivot_wider(names_from = Model, values_from = PValue)

# Print the pivot table
print(pivot_table)

library(fmsb)

# Prepare the data for radar chart
# radar_data <- as.data.frame(t(pivot_table[-1]))
# colnames(radar_data) <- pivot_table$lpvalue
# radar_data <- rbind(rep(max(radar_data), ncol(radar_data)), rep(min(radar_data), ncol(radar_data)), radar_data)
# 
# # Plot the radar chart
# radarchart(radar_data, axistype = 1,
#            title = "Radar Chart of Models by P-values",
#            pcol = rainbow(nrow(radar_data)-2), pfcol = rainbow(nrow(radar_data)-2, alpha=0.5), plwd = 2)
# legend(x = 1.1, y = 1.1, legend = rownames(radar_data)[-c(1,2)], col = rainbow(nrow(radar_data)-2), lty = 1, lwd = 2)
# 


library(ggplot2)
library(reshape2)
library(readxl)

pivot_table <- read_excel("C:/Users/vizama/Documents/1st paper/Real data Example/Canada weather reversed.xlsx")

# Melt the pivot table for plotting
pivot_long <- melt(pivot_table, id.vars = "lpvalue")

# Plot the facet plot
ggplot(pivot_long, aes(x = factor(lpvalue), y = value, fill = factor(lpvalue))) +
  geom_bar(stat = "identity") +
  facet_wrap(~ variable) +
  labs(title = "P-values by lp_value for Each Model",
       x = "lp_value", y = "PValue") +
  theme_minimal()




# Install ggpattern if you haven't already
# install.packages("remotes")
# remotes::install_github("coolbutuseless/ggpattern")

library(ggplot2)
library(reshape2)
library(ggpattern) 

# Create the bar plot with patterns for black-and-white printing

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 8),
  axis.title.y = element_text(size = 12))


ggplot(pivot_long, aes(x = factor(lpvalue,
                                  levels = c(1,2, 3, 4, 5, 10, 20,0))
                         , y = value, fill = variable, pattern = variable)) +
  geom_bar_pattern(stat = "identity", position = "dodge", color = "black",
                   pattern_density = 0.4,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02
  ) +
  scale_fill_grey(start = 0.1, end = 1) +  
  scale_pattern_manual(values = c("none", "crosshatch", "stripe", "circle", "none")) +  # Custom patterns
  labs(title = "",
       x = expression(paste("Size of p in"," ", L[italic(p)],"-", 'norm')),
       y = expression(paste(italic(p),'-',values)), fill = " ", pattern = " ") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.7, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE),
         pattern = guide_legend(nrow = 1, byrow = TRUE))+
  geom_hline(yintercept=0.05, linetype="dashed", color = "black",
             linewidth=1)+ My_Theme+ 
  scale_x_discrete(labels = c('0' = expression(infinity)))



ggplot(pivot_long, aes(x = factor(lpvalue), y = value, fill = variable, pattern = variable)) +
  geom_bar_pattern(stat = "identity", position = "dodge", color = "black",
                   pattern_density = 0.4,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02
  ) +
  scale_fill_grey(start = 0.1, end = 1) +  
  scale_pattern_manual(values = c("none", "crosshatch", "stripe", "circle", "none")) +  # Custom patterns
  labs(title = "",
       x = expression(paste("Size of p in"," ", L[italic(p)],"-", 'norm')),
       y = expression(paste(italic(p),'-',values)), fill = " ", pattern = " ") +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside =  c(.81, .95),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE),
         pattern = guide_legend(nrow = 1, byrow = TRUE))+
  geom_hline(yintercept=0.05, linetype="dashed", color = "black",
             linewidth=1)+ My_Theme



ggplot(pivot_long, aes(x = factor(lpvalue,
                                  levels = c(1,2, 3, 4, 5, 10, 20,0)),
                       y = value,
                       fill = variable, pattern = variable)) +
  geom_bar_pattern(stat = "identity", position = "dodge", color = "black",
                   pattern_density = 0.4,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02
  ) +
  scale_fill_grey(start = 0.1, end = 1) +  
  scale_pattern_manual(values = c("none", "crosshatch", "stripe", "circle", "none")) +  # Custom patterns
  labs(title = "",
       x = expression(paste("Size of p in"," ", L[italic(p)],"-", 'norm')),
       y = expression(paste(italic(p),'-',values)), fill = " ", pattern = " ") +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside =  c(.81, .95),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE),
         pattern = guide_legend(nrow = 1, byrow = TRUE))+
  geom_hline(yintercept=0.05, linetype="dashed", color = "black",
             linewidth=1)+ My_Theme + 
  scale_x_discrete






