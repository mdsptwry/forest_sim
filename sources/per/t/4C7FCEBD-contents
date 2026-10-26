# recruitment.R

# Method: batched k-d tree search (RANN::nn2) finds candidate parents for
# all gaps at once, then the Gumbel-max trick picks one weighted-random parent

recruit_trees <- function(forest, sim_config, species_params) {
  dead_idx  <- which(!forest$alive)
  n_dead    <- length(dead_idx)
  if (n_dead == 0) return(forest)

  alive_idx <- which(forest$alive)
  n_alive   <- length(alive_idx)
  if (n_alive == 0) return(forest)  # no possible parents left; leave gaps as-is

  fec_power       <- sim_config$recruitment$fecundity_power
  neighbor_radius <- sim_config$recruitment$neighbor_radius
  recruit_dbh     <- sim_config$recruitment$recruit_dbh

  global_crown_a  <- median(species_params$a,  na.rm = TRUE)
  global_crown_b  <- median(species_params$b,  na.rm = TRUE)
  global_height_a <- median(species_params$ta, na.rm = TRUE)
  global_height_b <- median(species_params$tb, na.rm = TRUE)
  global_height_k <- median(species_params$tk, na.rm = TRUE)

  params_lookup <- species_params |>
    select(species, crown_a = a, crown_b = b, height_a = ta, height_b = tb, height_k = tk)

  dead_coords  <- as.matrix(forest[dead_idx,  c("x", "y")])
  alive_coords <- as.matrix(forest[alive_idx, c("x", "y")])

  k_use <- min(200, n_alive)
  nn <- RANN::nn2(alive_coords, dead_coords, k = k_use)

  dead_rep  <- rep(seq_len(n_dead), times = k_use)
  nbr_local <- as.vector(nn$nn.idx)
  dist_vec  <- as.vector(nn$nn.dists)

  alive_dbh     <- forest$dbh_cm[alive_idx]
  alive_crown_r <- forest$crown_radius_m[alive_idx]

  nbr_dbh    <- alive_dbh[nbr_local]
  nbr_crownr <- alive_crown_r[nbr_local]
  sigma      <- pmax(nbr_crownr / 3.034, 1e-6)

  # Fecundity x dispersal-kernel weighting
  # formula from model: F = DBH^fecundity_power, k(d) = exp(-d^2 / 2*sigma^2).
  
  fec    <- nbr_dbh ^ fec_power
  kernel <- exp(-(dist_vec^2) / (2 * sigma^2))
  weight <- fec * kernel

  # if a gap has no candidates within radius, fall back to its nearest available candidates
  within_radius  <- dist_vec <= neighbor_radius
  has_any_within <- tapply(within_radius, dead_rep, any)
  keep <- within_radius | !has_any_within[dead_rep]

  dead_rep  <- dead_rep[keep]
  nbr_local <- nbr_local[keep]
  weight    <- weight[keep]

  # Gumbel-max trick: pick one weighted-random winner per gap
  weight[weight <= 0 | !is.finite(weight)] <- 1e-12
  score <- log(weight) - log(-log(runif(length(weight))))

  ord <- order(dead_rep, -score)
  sorted_dead_rep <- dead_rep[ord]
  winner <- ord[!duplicated(sorted_dead_rep)]

  parent_local  <- nbr_local[winner]
  parent_global <- alive_idx[parent_local]

  parent_species <- forest$species[parent_global]
  parent_x       <- forest$x[parent_global]
  parent_y       <- forest$y[parent_global]
  parent_crownr  <- forest$crown_radius_m[parent_global]

  # Recruit location: Gaussian dispersal centered on the parent, sigma
  # derived from the parent's own crown radius 
  sigma_dispersal <- pmax(parent_crownr / 3.034, 1e-6)
  x_new <- parent_x + rnorm(n_dead, 0, sigma_dispersal)
  y_new <- parent_y + rnorm(n_dead, 0, sigma_dispersal)

  # Boundary reflection — identical to the original.
  x_new <- ifelse(x_new < 0, abs(x_new), x_new)
  x_new <- ifelse(x_new > sim_config$plot_width,  2 * sim_config$plot_width  - x_new, x_new)
  y_new <- ifelse(y_new < 0, abs(y_new), y_new)
  y_new <- ifelse(y_new > sim_config$plot_height, 2 * sim_config$plot_height - y_new, y_new)

  # Species/allometry inheritance from the chosen parent, with the same
  # global-median fallback used elsewhere for species without a real fit 
  lookup_idx <- match(parent_species, params_lookup$species)
  new_crown_a  <- ifelse(is.na(lookup_idx), global_crown_a,  params_lookup$crown_a[lookup_idx])
  new_crown_b  <- ifelse(is.na(lookup_idx), global_crown_b,  params_lookup$crown_b[lookup_idx])
  new_height_a <- ifelse(is.na(lookup_idx), global_height_a, params_lookup$height_a[lookup_idx])
  new_height_b <- ifelse(is.na(lookup_idx), global_height_b, params_lookup$height_b[lookup_idx])
  new_height_k <- ifelse(is.na(lookup_idx), global_height_k, params_lookup$height_k[lookup_idx])

  next_id <- suppressWarnings(max(forest$id, na.rm = TRUE))
  if (!is.finite(next_id)) next_id <- 0
  new_ids <- seq(next_id + 1, next_id + n_dead)

  forest$id[dead_idx]       <- new_ids
  forest$species[dead_idx]  <- parent_species
  forest$x[dead_idx]        <- x_new
  forest$y[dead_idx]        <- y_new
  forest$dbh_cm[dead_idx]   <- recruit_dbh
  forest$age_days[dead_idx] <- 0
  forest$alive[dead_idx]    <- TRUE
  forest$crown_a[dead_idx]  <- new_crown_a
  forest$crown_b[dead_idx]  <- new_crown_b
  forest$height_a[dead_idx] <- new_height_a
  forest$height_b[dead_idx] <- new_height_b
  forest$height_k[dead_idx] <- new_height_k
  # New recruits haven't been through calculate_shading() yet 
  # this gives them a defined value of 0, they'll get a real value next time shading recomputes.
  forest$shading[dead_idx] <- 0
  forest$overlap[dead_idx] <- 0

  forest
}
