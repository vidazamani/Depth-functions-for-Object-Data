## Loading data
library(fda)

data(CanadianWeather)
x = CanadianWeather
CanadianWeather$dailyAv[,2,]


library(fda.usc)
###### Computing L_2 Distance 
Dismatrix = metric.lp(t(x$dailyAv[,,1]), lp = 2)




EMD <- function(D){
  
  n <- nrow(D)
  maximum = c()
  lens_depth = c()
  
 for (p in 1:n) {
   s = 0
   for (i in 1:n) {
     for (j in 1:n) {
       
       if (i<j) {
         maximum[i] = max(D[i,p],D[j,p])
         if (D[i,j]>maximum[i]) {
           s = s+1
         }
         
       }
       
     }
     
   }
   
   
   lens_depth[p]  = (1/choose(n,2))*s
   
 }
  
 return(lens_depth)
  
}


EMD(Dismatrix)




plot(EMD(Dismatrix),
     msd(Dismatrix), asp = 1)







