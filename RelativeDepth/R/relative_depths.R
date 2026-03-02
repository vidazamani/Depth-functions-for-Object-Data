################# Relative Depth Functions ############################
#' Modified Lens Depth (Relative Version)
#'
#' Computes the relative modified lens depth of observations
#' with respect to a reference subset.
#'
#' @param D A square distance matrix.
#' @param ref_idx Integer vector of reference indices.
#'
#' @return A numeric vector of depth values.
#' @export
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

#' Modified Metric Halp-space Depth (Relative Version)
#'
#' Computes the relative modified Metric Halp-space Depth of observations
#' with respect to a reference subset.
#'
#' @param D A square distance matrix.
#' @param ref_idx Integer vector of reference indices.
#'
#' @return A numeric vector of depth values.
#' @export
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


#' Modified Metric Spatial Depth (Relative Version)
#'
#' Computes the relative modified Metric Spatial Depth of observations
#' with respect to a reference subset.
#'
#' @param D A square distance matrix.
#' @param ref_idx Integer vector of reference indices.
#'
#' @return A numeric vector of depth values.
#' @export
MSD_relative <- function(D, ref_idx) {

  n_total <- nrow(D)
  n_ref   <- length(ref_idx)

  res <- rep(0, n_total)

  for (k in 1:n_total) {

    res_now <- 0

    for (i in ref_idx) {
      for (j in ref_idx) {

        if (D[k,i] > 1e-6 && D[k,j] > 1e-6) {

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



#' Modified Metric Oja Depth 2D (Relative Version)
#'
#' Computes the relative modified Metric Oja Depth (2D) of observations
#' with respect to a reference subset.
#'
#' @param D A square distance matrix.
#' @param ref_idx Integer vector of reference indices.
#' @return A numeric vector of depth values.
#' @export
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


#' Modified Metric Oja Depth 3D (Relative Version)
#'
#' Computes the relative modified Metric Oja Depth (3D) of observations
#' with respect to a reference subset.
#'
#' @param D A square distance matrix.
#' @param ref_idx Integer vector of reference indices.
#' @return A numeric vector of depth values.
#' @export
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

          area <- area + sqrt(detS + 4*D[i,w]^2*D[j,w]^2*D[k,w]^2)
        }
      }
    }

    ojadepth[w] <- 1 / (1 + (1/(n_ref*(n_ref-1)*(n_ref - 2)/6)) * area)
  }

  return(ojadepth)
}
