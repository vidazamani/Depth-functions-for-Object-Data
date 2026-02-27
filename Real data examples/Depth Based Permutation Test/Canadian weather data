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


my_col <- brewer.pal(4, "Dark2")
plot(ft_s, col = my_col[factor(region)], lty = 1, xlab = "Time (days)", ylab = "Average temperature")
legend(1, 20, legend = c("Arctic", "Atlantic", "Continental", "Pacific"), lwd = 1, col = my_col, cex = 0.8)


my_col_2 <- my_col
my_col_2[4] <- my_col[2]
my_col_2[1] <- my_col[3]

plot(ft_s, col = my_col_2[factor(region)], lty = 1, xlab = "Time (days)", ylab = "Average temperature")
legend(1, 20, legend = c("Arctic or Continental", "Atlantic or Pacific"), 
       lwd = 1,
       col = my_col_2[c(1, 2, 3)],
       cex = 0.8)


# The mean function corresponding to a sample of curves stored in an fd-object can be computed with mean.fd
ft_mean <- mean.fd(ft_s)
plot(ft_s, col = 1, lty = 2, xlab = "Time (days)", ylab = "Average temperature")
plot(ft_mean, add = TRUE, col = 2, lwd = 3, lty = 1)



#Divide the stations into two groups ("Atlantic and Pacific" vs. "Continental and Arctic") and compute the mean and +-2 SD curves separately for the two groups.

g1 <- region %in% c("Atlantic", "Pacific")
g2 <- region %in% c("Continental", "Arctic")


ftg1 <- ft_s[g1]
ftg1_mean <- mean.fd(ftg1)
ftg1_sd <- sd.fd(ftg1)


ftg2 <- ft_s[g2]
ftg2_mean <- mean.fd(ftg2)
ftg2_sd <- sd.fd(ftg2)



plot(ftg1_mean, ylim = c(-25, 20), xlab = "Time (days)", ylab = "Mean temperature across stations")
plot(ftg2_mean, col = 2, add = TRUE)
legend(1, 20, legend = c("Atlantic and Pacific", "Continental and Arctic"), lwd = 1, col = c(1, 2), cex = 0.8)


plot(ftg1_sd, ylim = c(0, 10), xlab = "Time (days)", ylab = "SD of temperature across stations")
plot(ftg2_sd, col = 2, add = TRUE)
legend(1, 10, legend = c("Atlantic and Pacific", "Continental and Arctic"), lwd = 1, col = c(1, 2), cex = 0.8)





#################################################################################################
######### Data Preparation

# ap <- CanadianWeather$dailyAv[,,1][,g1]

# ca <- CanadianWeather$dailyAv[,,1][,g2]


###### Computing L_2 Distance which is needed to be used in depth functions to find the deepest object
# Dismatrixap = metric.lp(t(ap), lp = 2)
Dismatrixap2 <- metric.lp(fdata(ftg1), lp = 2)


# Dismatrixca = metric.lp(t(ca), lp = 2)
Dismatrixca2 <- metric.lp(fdata(ftg2), lp = 2)


#Difference in means for real observed data

ordap <- order(MOD3(Dismatrixap2),decreasing = TRUE)[1]
deepestap <- fdata(ftg1)[ordap]

ordca <- order(MOD3(Dismatrixca2),decreasing = TRUE)[1]
deepestca <- fdata(ftg2)[ordca]




# clean up
par(mar = c(5.1,4.1,5,5))


### Smoothed Data

plot(ft_s, col = my_col[factor(g1)], lty = 1, xlab = "Time (days)", ylab = "Average temperature")

plot(ft_s[(1:35)[g1][ordap]], add = TRUE,
     col = "black", lty = 1,lwd = 3,
     xlab = "Time (days)", ylab = "Average temperature")
plot(ft_s[(1:35)[g2][ordca]], col = "Blue", lwd = 3,add = TRUE,
     lty = 1, xlab = "Time (days)", ylab = "Average temperature")



#### Raw Data

# matplot(day.5, CanadianWeather$dailyAv[, , "Temperature.C"],
#         type="l", axes=FALSE, xlab="", ylab="Mean Temperature (deg C)",
#         col = my_col[factor(g1)])
# axis(2, las=1)

highlight_lines <- c((1:35)[g1][ordap], (1:35)[g2][ordca])
# 
# colors <- c("blue", "green") # Colors for highlighted lines
# line_types <- c(2, 4)      # Line types for highlighted lines (dashed and dotted)
# line_widths <- c(2, 3)     # Line widths for highlighted lines

# Use matlines to add the highlighted lines
# matlines(day.5, CanadianWeather$dailyAv[, , "Temperature.C"][, highlight_lines], 
#          col = colors,
#          lty = line_types,
#          lwd = line_widths)


# Label the horizontal axis with the month names
# axis(1, monthBegin.5, labels=FALSE)
# axis(1, monthEnd.5, labels=FALSE)
# axis(1, monthMid, monthLetters, tick=FALSE)

# Add the names of the weather stations
p  = c(colnames(CanadianWeather$dailyAv[, , "Temperature.C"][, highlight_lines]))
mtext(p, side=4,
      at=CanadianWeather$dailyAv[365, p, "Temperature.C"],
      las=1)



pvalue_canada_weather <- function(lp_value, iter, depth_fun) {
  
  distribution=c()
  
  for (i in 1:iter) {
    
    fake_labels <- sample(c(rep("ap",20),rep("ca",15)), 35, replace = FALSE)
    
    ap_fake <- ft_s[fake_labels %in% c("ap")]
    
    ca_fake <- ft_s[fake_labels %in% c("ca")]
    
    
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
  
  Dismatrixap_real <- metric.lp(fdata(ftg1), lp = lp_value)
  
  Dismatrixca_real <- metric.lp(fdata(ftg2), lp = lp_value)
  
  
  ordap_r <- order(depth_fun(Dismatrixap_real),decreasing = TRUE)[1]
  deepestap_real <- fdata(ftg1)[ordap_r]
  
  ordca_r <- order(depth_fun(Dismatrixca_real),decreasing = TRUE)[1]
  deepestca_real <- fdata(ftg2)[ordca_r]
  
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
  data.frame(Model = combinations$ModelName[i], lpvalue = lp_val, PValue = pvalue_canada_weather(lp_val, iter, model))
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

# # Prepare the data for radar chart
# radar_data <- as.data.frame(t(pivot_table[-1]))
# colnames(radar_data) <- pivot_table$lpvalue
# radar_data <- rbind(rep(max(radar_data), ncol(radar_data)), rep(min(radar_data), ncol(radar_data)), radar_data)
# 
# # Plot the radar chart
# radarchart(radar_data, axistype = 1,
#            title = "Radar Chart of Models by P-values",
#            pcol = rainbow(nrow(radar_data)-2), pfcol = rainbow(nrow(radar_data)-2, alpha=0.5), plwd = 2)
# legend(x = 1.1, y = 1.1, legend = rownames(radar_data)[-c(1,2)], col = rainbow(nrow(radar_data)-2), lty = 1, lwd = 2)



library(ggplot2)
library(reshape2)
library(readxl)

pivot_table <- read_excel("C:/Users/vizama/Documents/1st paper/Real data Example/canada weather.xlsx")

# Melt the pivot table for plotting
pivot_long <- melt(pivot_table, id.vars = "lpvalue")

# Plot the facet plot
ggplot(pivot_long, aes(x = factor(lpvalue), y = value, fill = factor(lpvalue))) +
  geom_bar(stat = "identity") +
  facet_wrap(~ variable) +
  labs(title = "P-values by lp_value for Each Model",
       x = "lp_value", y = "PValue") +
  theme_minimal()


# Plot the bar plot
ggplot(pivot_long, aes(x = factor(lpvalue), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "P-values by lp_value and Model",
       x = "lp_value", y = "P-values", fill = "Model") +
  theme_minimal()



# Install ggpattern if you haven't already
# install.packages("remotes")
# remotes::install_github("coolbutuseless/ggpattern")

library(ggplot2)
library(reshape2)
library(ggpattern) 

# Create the bar plot with patterns for black-and-white printing
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
       x = expression(paste("Size of p in"," ", L[p],"-", 'norm')),
       y = "P-value", fill = " ", pattern = " ") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE),
         pattern = guide_legend(nrow = 1, byrow = TRUE))


My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 8),
  axis.title.y = element_text(size = 12))

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
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.7, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE),
         pattern = guide_legend(nrow = 1, byrow = TRUE))+
  geom_hline(yintercept=0.05, linetype="dashed", color = "black",
             linewidth=1)+ My_Theme



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
  theme(legend.position = "bottom",
        legend.key.size = unit(0.7, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE),
         pattern = guide_legend(nrow = 1, byrow = TRUE))+
  geom_hline(yintercept=0.05, linetype="dashed", color = "black",
             linewidth=1)+ My_Theme


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





# Plot the line plot
# ggplot(pivot_long, aes(x = lpvalue, y = value, color = variable, group = variable)) +
#   geom_line() +
#   geom_point() +
#   labs(title = "P-values Across lp_values for Each Model",
#        x = "lp_value", y = "P-values", color = "Model") +
#   theme_minimal()
# 
# 
# # Plot the heatmap
# ggplot(pivot_long, aes(x = variable, y = factor(lpvalue), fill = value)) +
#   geom_tile(color = "white") +
#   scale_fill_gradient(low = "white", high = "blue") +
#   labs(title = "Heatmap of P-values by Model and lp_value",
#        x = "Model", y = "lp_value", fill = "P-values") +
#   theme_minimal()
# 
# 
# # Plot the black-and-white heatmap
# ggplot(pivot_long, aes(x = variable, y = factor(lpvalue), fill = value)) +
#   geom_tile(color = "black") +  # Add black borders to the tiles
#   scale_fill_gradient(low = "white", high = "black", na.value = "gray90") +  # Shades of gray
#   geom_text(aes(label = round(value, 2)), color = "black", size = 3) +  # Add text labels
#   labs(title = "Heatmap of PValues by Model and lp_value (B&W)",
#        x = "Model", y = "lp_value", fill = "PValue") +
#   theme_minimal(base_size = 14) +  # Increase base text size for visibility
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
#     panel.grid = element_blank(),  # Remove gridlines
#     plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
#   )
# 
# ggplot(pivot_long, aes(x = variable, y = factor(lpvalue), fill = value)) +
#   geom_tile(color = "black") +  # Add black borders to tiles
#   scale_fill_gradient(low = "white", high = "black", na.value = "gray90") +  # Shades of gray
#   geom_text(aes(label = round(value, 2)), color = "black", size = 3) +  # Add text labels
#   geom_text(aes(label = ifelse(value > 0.05, "*", "")), vjust = -1.5, color = "red", size = 3) +  # Annotate with asterisks
#   labs(title = "Heatmap of PValues by Model and lp_value with Annotations",
#        x = "Model", y = "lp_value", fill = "PValue") +
#   theme_minimal(base_size = 14) +  # Increase base text size for visibility
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
#     panel.grid = element_blank(),  # Remove gridlines
#     plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
#   )
# 
# ggplot(pivot_long, aes(x = variable, y = factor(lpvalue), fill = value)) +
#   geom_tile(color = "black") +
#   scale_fill_gradient(low = "white", high = "black", na.value = "gray90") +
#   geom_text(aes(label = round(value, 2)), color = "black", size = 4, fontface = "bold") +  # Bold labels for clarity
#   labs(title = "Enhanced Heatmap of PValues", x = "Model", y = "lp_value", fill = "PValue") +
#   theme_minimal(base_size = 16) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "italic"),  # Italicize x-axis labels
#     axis.text.y = element_text(size = 12, face = "italic"),  # Italicize y-axis labels
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  # Larger, centered title
#     panel.grid = element_blank()
#   )
# 
# ggplot(pivot_long, aes(x = variable, y = factor(lpvalue), fill = value)) +
#   geom_tile(color = "black") +  # Add black borders to tiles
#   scale_fill_gradient(low = "white", high = "black", limits = c(0, 1), na.value = "gray90") +  # Normalize color scale from 0 to 1
#   geom_text(aes(label = round(value, 2)), color = "black", size = 3) +  # Add text labels
#   labs(title = "Heatmap of PValues by Model and lp_value (Normalized to 0-1)",
#        x = "Model", y = "lp_value", fill = "PValue") +
#   theme_minimal(base_size = 14) +  # Increase base text size for visibility
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
#     panel.grid = element_blank(),  # Remove gridlines
#     plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
#   )
# 
# # Plot the heatmap with full PValue range from 0 to 1
# ggplot(pivot_long, aes(x = variable, y = factor(lpvalue), fill = value)) +
#   geom_tile(color = "white") +
#   scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1), na.value = "grey50") +
#   labs(title = "Heatmap of PValues by Model and lp_value",
#        x = "Model", y = "lp_value", fill = "PValue") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # Convert Model to a factor with specific levels for control over order
# pivot_long$variable <- factor(pivot_long$variable, levels = unique(pivot_long$variable))
# 
# # Plot the bar plot optimized for black-and-white printing
# ggplot(pivot_long, aes(x = factor(lpvalue), y = value, fill = variable, pattern = variable)) +
#   geom_bar(stat = "identity", position = "dodge", color = "black") +
#   scale_fill_grey(start = 0.3, end = 0.9) +
#   labs(title = "PValues by lp_value and Model",
#        x = "lp_value", y = "PValue", fill = "Model") +
#   theme_minimal() +
#   theme(legend.position = "bottom",
#         legend.key.size = unit(0.5, "cm"),
#         axis.text.x = element_text(angle = 45, hjust = 1)) +
#   guides(fill = guide_legend(nrow = 2, byrow = TRUE))




hist(distribution, breaks=10, col='grey', main="Permutation Distribution", las=1, xlab='')
abline(v=observed, lwd=3, col="red")
