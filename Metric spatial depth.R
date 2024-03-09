#----------------

## Loading data
library(fda)
data(CanadianWeather)
x = CanadianWeather
CanadianWeather$dailyAv[1,,]

## 1.  Plot (latitude & longitude) of stations by region
with(CanadianWeather, plot(-coordinates[, 2], coordinates[, 1], type='n',
                           xlab="West Longitude", ylab="North Latitude",
                           axes=FALSE) )
Wlon <- pretty(CanadianWeather$coordinates[, 2])
axis(1, -Wlon, Wlon)
axis(2)

rgns <- 1:4
names(rgns) <- c('Arctic', 'Atlantic', 'Continental', 'Pacific')
Rgns <- rgns[CanadianWeather$region]
with(CanadianWeather, points(-coordinates[, 2], coordinates[, 1],
                             col=Rgns, pch=Rgns) )
legend('topright', legend=names(rgns), col=rgns, pch=rgns)



##
## 2.  Plot dailyAv[, 'Temperature.C'] for 4 stations
##
data(CanadianWeather)
# Expand the left margin to allow space for place names
op <- par(mar=c(5, 4, 4, 5)+.1)
# Plot
stations <- c("Pr. Rupert", "Montreal", "Edmonton", "Resolute")
matplot(day.5, CanadianWeather$dailyAv[, stations, "Temperature.C"],
        type="l", axes=FALSE, xlab="", ylab="Mean Temperature (deg C)")
axis(2, las=1)
# Label the horizontal axis with the month names
axis(1, monthBegin.5, labels=FALSE)
axis(1, monthEnd.5, labels=FALSE)
axis(1, monthMid, monthLetters, tick=FALSE)
# Add the monthly averages
matpoints(monthMid, CanadianWeather$monthlyTemp[, stations])
# Add the names of the weather stations
mtext(stations, side=4,
      at=CanadianWeather$dailyAv[365, stations, "Temperature.C"],
      las=1)
# clean up
par(op)



# Metric spatial depth

# Below, D contains non-squared distances

# Main function, the metric spatial depth of a sample

msd <- function(D){
  
  n <- nrow(D)
  
  res <- rep(0, n)
  
  for(k in 1:n){
    
    res_now <- 0
    
    d_row <- D[k, ]
    
    index_set <- setdiff(1:n, k)
    
    index_set <- setdiff(index_set, which(d_row < 1e-6))
    
    for(i in index_set){
      
      for(j in index_set){
        
        temp <- d_row[i]/d_row[j]
        
        res_now <- res_now + temp + 1/temp - D[i, j]^2/(d_row[i]*d_row[j])
        
      }
      
    }
    
    res[k] <- res_now
    
  }
  
  res <- res/n^2
  
  1 - 0.5*res
  
}

# The depth of a new point

# D = n x n distance matrix of a sample

# d = n-vector of distances between a point and the sample

msd_test <- function(D, d){
  
  n <- nrow(D)
  
  res <- 0
  
  d_row <- d
  
  index_set <- setdiff(1:n, which(d_row < 1e-6))
  
  for(i in index_set){
    
    for(j in index_set){
      
      temp <- d_row[i]/d_row[j]
      
      res <- res + temp + 1/temp - D[i, j]^2/(d_row[i]*d_row[j])
      
    }
    
  }
  
  res <- res/n^2
  
  1 - 0.5*res
  
}

# Kernel depth of a new point

# (the code assumes that k(x, x) = 1!)

# K = n x n kernel matrix

# k = n-vector of kernel values between a new point and a sample

msd_kernel_test <- function(K, k){
  
  n <- nrow(K)
  
  res <- 0
  
  d_row <- diag(K) - 2*k + 1
  
  index_set <- setdiff(1:n, which(d_row < 1e-6))
  
  d_row <- sqrt(d_row)
  
  for(i in index_set){
    
    for(j in index_set){
      
      res <- res + (K[i, j] - k[i] - k[j] + 1)/(d_row[i]*d_row[j])
      
    }
    
  }
  
  res <- res/n^2
  
  1 - res
  
}

# ISOMAP-style non-linear depth of a collection of test points

library(vegan)

# D = distance matrix of the training set

# D0 = distances between the training set and the test set (n x n0)

# k = neighbours

msd_isomap_test <- function(D, D0, k){
  
  n <- ncol(D)
  
  n0 <- ncol(D0)
  
  # Remove the non-neighbour distances
  
  d_max <- max(cbind(D, D0))
  
  D_temp <- as.matrix(D)
  
  diag(D_temp) <- NA
  
  is.na(D_temp) <- apply(D_temp, 2, function(x) x > x[order(x, na.last=TRUE)[k]])
  
  D_graph <- pmax(as.dist(D_temp), as.dist(t(D_temp)), na.rm = TRUE)
  
  res_full <- rep(0, n0)
  
  for(i in 1:n0){
    
    # Graft the test point and obtain distances
    
    d0 <- D0[, i]
    
    is.na(d0) <- (d0 > d0[order(d0, na.last = TRUE)[k]])
    
    D_aug <- as.dist(rbind(cbind(as.matrix(D_graph), d0), c(d0, 0)))
    
    D_aug[is.na(D_aug)] <- d_max + 2
    
    D_aug_graph <- stepacross(D_aug, toolong = d_max + 1)
    
    D_full <- as.matrix(D_aug_graph)[1:n, 1:n]
    
    d_row <- as.matrix(D_aug_graph)[-(n + 1), (n + 1)]
    
    index_set <- setdiff(1:n, which(d_row < 1e-6))
    
    res <- 0
    
    for(ell in index_set){
      
      for(j in index_set){
        
        temp <- d_row[ell]/d_row[j]
        
        res <- res + temp + 1/temp - D_full[ell, j]^2/(d_row[ell]*d_row[j])
        
      }
      
    }
    
    res_full[i] <- res
    
  }
  
  res_full <- res_full/n^2
  
  1 - 0.5*res_full
  
}