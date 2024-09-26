#### Real example with Canadian weather data ####


##### All needed source and libraries ########
library(fda)
library(fda.usc) ### For computing distance
library(RColorBrewer)


source("C:/1stPaper/Codes/All Metric depth functions.R")


###################################################################################################
########## In this part we just want to show why this data can be divided to two groups ###########
# The full temperature data set from the Canadian weather data example
y <- CanadianWeather$dailyAv[, , 1]
my_basis <- create.fourier.basis(c(0, 365), 100)
ft <- Data2fd(1:365, y, my_basis)
ft_s <- smooth.fd(ft, fdPar(ft, 2, 1e4))


# Each of the weather stations belongs to one of four regions:
# Atlantic, Continental, Pacific, Arctic

region <- CanadianWeather$region


my_col <- brewer.pal(4, "Dark2")

my_col_2 <- my_col
my_col_2[4] <- my_col[2]
my_col_2[1] <- my_col[3]



par(mfrow = c(1,3))

# Define line types for each region
line_types <- c(1, 2, 1, 2)  # You can choose different types of lines: 1 = solid, 2 = dashed, 3 = dotted
line_types[4] <- line_types[2]
line_types[1] <- line_types[3]



# Plot with color and line types
plot(ft_s, col = my_col_2[factor(region)],
     lty = line_types[factor(region)], 
     lwd = 1.5,
     xlab = "Time (days)", ylab = "Daily Average temperature",
     cex.lab=0.75)

# Add legend with color and line types
legend(-20, 24, legend = c("Atlantic or Pacific", "Arctic or Continental"), 
       lwd = 1.5, col = my_col_2[c(2, 1, 3)], 
       lty = line_types[c(2, 1, 3)], 
       cex = 0.5,
       bty = 'n')






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



plot(ft_s, col = 'white', lty = 1, xlab = "Time (days)",
     ylab = "the deepest temperature of each group", cex.lab = 0.75)

plot(ft_s[(1:35)[g1][ordap]], add = TRUE,
     col = 2, lty = 2,lwd = 3,
     xlab = "Time (days)")
plot(ft_s[(1:35)[g2][ordca]], col = 1, lwd = 1.5,add = TRUE,
     lty = 1, xlab = "Time (days)")



legend(-20, 22, legend = c("Atlantic or Pacific", "Arctic or Continental"),
       lwd = c(3,1.5),
       col = c(2, 1),
       cex = 0.5,
       lty = c(2,1),
       bty = "n")





#Divide the stations into two groups ("Atlantic and Pacific" vs. "Continental and Arctic") and compute the mean and +-2 SD curves separately for the two groups.

g1 <- region %in% c("Atlantic", "Pacific")
g2 <- region %in% c("Continental", "Arctic")


ftg1 <- ft_s[g1]
ftg1_mean <- mean.fd(ftg1)


ftg2 <- ft_s[g2]
ftg2_mean <- mean.fd(ftg2)



plot(ftg1_mean, ylim = c(-25, 20), xlab = "Time (days)",
     ylab = "mean and the deepest temperature of each group", lwd = 3,
     lty = 4, col = 3,cex.lab=0.75)
plot(ftg2_mean, col = 4, add = TRUE, lwd = 2,
     lty = 1)



plot(ft_s[(1:35)[g1][ordap]], add = TRUE,
     col = 2, lty = 3,lwd = 4,
     xlab = "Time (days)")
plot(ft_s[(1:35)[g2][ordca]], col = 1, lwd = 2,add = TRUE,
     lty = 2, xlab = "Time (days)")




legend(-20, 22, legend = c("AP mean", "AC mean",
                           'AP deepest', 'AC deepest'),
       lwd = c(3, 2, 4, 2),
       col = c(3,4,2,1),
       cex = 0.5,
       lty = c(4,1,3,2),
       bty = "n")


########### Just Mean between each group


plot(ftg1_mean, ylim = c(-25, 20), xlab = "Time (days)",
     ylab = "mean temperature of each group", lwd = 3,
     lty = 4, col = 2,cex.lab=0.75)
plot(ftg2_mean, col = 1, add = TRUE, lwd = 2,
     lty = 1)


legend(-20, 22, legend = c('Atlantic or Pacific', 'Arctic or Continental'),
       lwd = c(3, 2),
       col = c(2,1),
       cex = 0.5,
       lty = c(4,1),
       bty = "n")

