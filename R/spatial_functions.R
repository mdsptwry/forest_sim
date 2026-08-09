# spatial_functions.R

# spatial operations for overlap and shading

# KNN search for faster overlap calculation

calculate_overlaps_knn <- function(data, k_neighbors = 50) {
  
  n <- nrow(data)
  
  if (n == 0) {
    return(data |>
             mutate(overlap = 0, area = pi * size^2, shading = 0))
  }
  
  if (n == 1) {
    return(data |>
             mutate(overlap = 0, area = pi * size^2, shading = 0))
  }
  
  # Step 1: Candidate pair search via k-d tree (RANN::nn2)
  # We search for each tree's k nearest neighbors. 
  # k_neighbors set to 50 for headroom (can be adjusted later)
  # If a tree has more true overlaps than k_neighbors, increase k_neighbors.
  k_use <- min(k_neighbors + 1, n)  # +1 because nn2 includes self as neighbor 0
  
  coords <- as.matrix(data[, c("x", "y")])
  nn_result <- RANN::nn2(coords, coords, k = k_use)
  
  # nn_result$nn.idx: n x k_use matrix of neighbor row-indices (col 1 = self)
  # nn_result$nn.dists: matching distances
  
  idx_self <- rep(seq_len(n), times = k_use - 1)
  idx_nbr  <- as.vector(nn_result$nn.idx[, 2:k_use])   # drop self column
  d_vec    <- as.vector(nn_result$nn.dists[, 2:k_use])
  
  # Step 2: Filter to pairs within plausible overlap range (sum of radii)
  r_self <- data$size[idx_self]
  r_nbr  <- data$size[idx_nbr]
  r_sum  <- r_self + r_nbr
  
  keep <- d_vec <= r_sum & d_vec > 0   # exclude zero-distance duplicates/self
  # drop pairs where the candidate radius exceeded k_neighbors search radius
  # (warn if k_neighbors is too small)
  idx_self <- idx_self[keep]
  idx_nbr  <- idx_nbr[keep]
  d_vec    <- d_vec[keep]
  
  # de-duplicate symmetric pairs (i,j) and (j,i) -> keep one direction only,
  # we'll still compute both g_i and g_j from that single row
  pair_key <- ifelse(idx_self < idx_nbr,
                     paste(idx_self, idx_nbr, sep = "_"),
                     paste(idx_nbr, idx_self, sep = "_"))
  dup <- duplicated(pair_key)
  idx_self <- idx_self[!dup]
  idx_nbr  <- idx_nbr[!dup]
  d_vec    <- d_vec[!dup]
  
  if (length(idx_self) == 0) {
    return(data |>
             mutate(overlap = 0, area = pi * size^2, shading = 0))
  }
  
  r1 <- data$size[idx_self]
  r2 <- data$size[idx_nbr]
  h1 <- data$height_m[idx_self]
  h2 <- data$height_m[idx_nbr]
  d  <- d_vec
  
  # Step 3: Fully vectorized circle-circle overlap area
  # Same three geometric cases as before, computed across all pairs at once.
  overlap <- numeric(length(d))
  
  case_contained <- d <= abs(r1 - r2)
  case_separate  <- d >= (r1 + r2)
  case_partial   <- !case_contained & !case_separate
  
  # case 1: one circle fully contains the other
  overlap[case_contained] <- pi * pmin(r1[case_contained], r2[case_contained])^2
  
  # case 2: no overlap
  overlap[case_separate] <- 0
  
  # case 3: partial lens-shaped overlap (vectorized lens formula)
  if (any(case_partial)) {
    r1p <- r1[case_partial]; r2p <- r2[case_partial]; dp <- d[case_partial]
    
    alpha <- 2 * acos(pmin(pmax((r1p^2 + dp^2 - r2p^2) / (2 * r1p * dp), -1), 1))
    beta  <- 2 * acos(pmin(pmax((r2p^2 + dp^2 - r1p^2) / (2 * r2p * dp), -1), 1))
    
    area1 <- 0.5 * r1p^2 * (alpha - sin(alpha))
    area2 <- 0.5 * r2p^2 * (beta  - sin(beta))
    
    overlap[case_partial] <- area1 + area2
  }
  
  # Step 4: Height-weighted asymmetric contributions (vectorized)
  g1 <- h2 / (h1 + h2)   # credited to idx_self (tree 1 / "i")
  g2 <- h1 / (h1 + h2)   # credited to idx_nbr  (tree 2 / "j")
  
  shade_contrib_1 <- overlap * g1
  shade_contrib_2 <- overlap * g2
  
  # Step 5: Aggregate per tree (vectorized via tapply, no group_by needed)
  contrib_ids   <- c(idx_self, idx_nbr)
  contrib_vals  <- c(shade_contrib_1, shade_contrib_2)
  
  agg <- tapply(contrib_vals, contrib_ids, sum)
  
  weighted_overlap <- numeric(n)
  weighted_overlap[as.integer(names(agg))] <- as.numeric(agg)
  
  # Step 6: Final shading index
  result <- data |>
    mutate(
      area    = pi * size^2,
      overlap = weighted_overlap,
      shading = pmin(weighted_overlap / area, 1)
    )
  
  return(result)
}
