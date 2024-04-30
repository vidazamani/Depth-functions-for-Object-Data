######### Metric Depth Functions  ######### 

# Lens depth for Object data
MLD <- function(D){
  
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



# Tukey's depth for object data (Metric Halfspace depth)
MHD <- function(D){
  
  n <- nrow(D)
  p = matrix(0,n,n)
  tu = c()
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      s = 0
      for (k in 1:n) {
        if (D[k,i]<= D[k,j]+1e-6){
          s=s+1
        }
      }
      
      p[i,j] = (1/n)*s
      
    }
  }
  
  for (y in 1:n) {
    Q = NULL
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        
        if (D[y,i]<=D[y,j]){
          Q = c(Q,p[i,j])
        }
      }
    }
    tu[y] = min(Q)
  }
  
  
  
  return(tu)
  
}

# Metric Spatial depth
MSD <- function(D){
  
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


# Metric Oja depth in 2 dimensional 
MOD2 <- function(D){
  
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


# Metric Oja depth in three-dimensional
MOD3 <- function(D){
  
  n <- nrow(D)
  B = array(0,c(n,n,n))
  ojadepth = c()
  area = c(0)
  
  
  for (w in 1:n) {
    
    for (i in 1:(n-1)) {
      for (j in i:n) {
        for (k in 1:n) {
          
          S = matrix(c((D[i,w])^2,
                       -0.5*((D[j,i])^2 - (D[j,w]^2) - (D[i,w])^2),
                       -0.5*((D[k,i])^2 - (D[k,w]^2) - (D[i,w])^2),
                       -0.5*((D[i,j])^2 - (D[i,w]^2) - (D[j,w])^2),
                       (D[j,w])^2,
                       -0.5*((D[k,j])^2 - (D[k,w]^2) - (D[j,w])^2),
                       -0.5*((D[i,k])^2 - (D[i,w]^2) - (D[k,w])^2),
                       -0.5*((D[j,k])^2 - (D[j,w]^2) - (D[k,w])^2),
                       (D[k,w])^2 )
                     ,3,3)
          
          B[i,j,k] = det(S)
          
        }
      }
    }
    B[B<0] = 0
    area[w] = sum(sqrt(B))
    ojadepth[w] = 1 / (1 + (1/(0.5*(n-1)*(n^2))*(area[w])))
  }
  
  return(ojadepth)
  
}