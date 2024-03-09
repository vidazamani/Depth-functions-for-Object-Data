library(ICtest)
library(abind)
library(CovTools)
library(ggplot2)
library(dplyr)


## Loading data
library(fda)

data(CanadianWeather)
x = CanadianWeather
CanadianWeather$dailyAv[,2,]


library(fda.usc)
###### Computing L_2 Distance 
## lp = 0 ----> infinity   lp = 2 ------> L_2 distance
Dismatrix = metric.lp(t(x$dailyAv[,,1]), lp = 2)













set.seed(12)
# Oja depth
Oja_v2 <- function(D){
  
  n <- nrow(D)
  p = matrix(0,n,n)
  ojadepth = c()
  area = c(0)
  
  
  for (k in 1:n) {
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        S = matrix(c((D[i,k])^2,
                     -0.5*((D[j,i])^2 - (D[j,k]^2) - (D[i,k])^2),
                     -0.5*((D[i,j])^2 - (D[i,k]^2) - (D[j,k])^2),
                     (D[j,k])^2)
                   ,2,2)
        
        p[i,j] = det(S)
        
      }
    }
    #area[k] = sum(upper.tri(p)*p)
    area[k] = sum(sqrt(p))
    ojadepth[k] = 1 / (1 + (1/(0.5*(n^2-n))*(area[k])))
  }
  
  return(ojadepth)
  
}


Oja_v2(Dismatrix)*1000


order(Oja_v2(Dismatrix)*1000)



