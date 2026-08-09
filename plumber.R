# plumber.R
# ============================================================

# Run it with:
#   plumber::plumb("plumber.R") |> plumber::pr_run(port = 8000)

library(plumber)
library(tidyverse)

source("R/allometry.R")
source("R/spatial_functions.R")
source("R/shading.R")
source("R/growth.R")
source("R/mortality.R")
source("R/recruitment.R")

# ---- helper: load a real BCI census file --------------------------------
load_census <- function(census_number) {
  path <- sprintf("data/bci.tree%d.rdata", census_number)
  if (!file.exists(path)) {
    stop(sprintf("Census file not found: %s", path))
  }
  e <- new.env()
  load(path, envir = e)
  obj_name <- sprintf("bci.tree%d", census_number)
  if (!exists(obj_name, envir = e)) {
    loaded_names <- ls(e)
    if (length(loaded_names) == 1) {
      obj_name <- loaded_names[1]
    } else {
      stop(sprintf(
        "Could not find object '%s' in %s. Objects found instead: %s",
        obj_name, path, paste(loaded_names, collapse = ", ")
      ))
    }
  }
  get(obj_name, envir = e)
}

# ---- Serve the front-end -----------------------------------------------

#* @assets ./www /
list()

# ---- Species allometry parameters (real BCI data) -------------------------
species_params <- read.csv("data/species_allometry_params.csv")

# ---- Full species pool (all ~307 species, from real census data) ----------
full_species_pool <- tryCatch({
  raw <- load_census(1)
  alive_rows <- raw[raw$status %in% c("A", "AD", "AM", "AR"), ]
  codes <- unique(alive_rows$sp)
  codes <- codes[!is.na(codes) & codes != ""]
  sort(codes)
}, error = function(e) {
  message("Could not load census 1 to build the full species pool (", conditionMessage(e),
          "). Falling back to the 162 species in species_allometry_params.csv.")
  sort(unique(species_params$species))
})
message(sprintf("Full species pool loaded: %d species (%d have real fitted allometry).",
                 length(full_species_pool), sum(full_species_pool %in% species_params$species)))

# ---- helper: resolve a set of species names into a full params table ------

resolve_species_subset <- function(names) {
  global_a  <- median(species_params$a,  na.rm = TRUE)
  global_b  <- median(species_params$b,  na.rm = TRUE)
  global_ta <- median(species_params$ta, na.rm = TRUE)
  global_tb <- median(species_params$tb, na.rm = TRUE)
  global_tk <- median(species_params$tk, na.rm = TRUE)

  tibble(species = names) |>
    left_join(species_params, by = "species") |>
    mutate(
      a  = ifelse(is.na(a),  global_a,  a),
      b  = ifelse(is.na(b),  global_b,  b),
      ta = ifelse(is.na(ta), global_ta, ta),
      tb = ifelse(is.na(tb), global_tb, tb),
      tk = ifelse(is.na(tk), global_tk, tk)
    )
}


# ---- helper: procedural initial layout ------------------------------------
# Grid-ish placement with soft cluster pull, so the demo is self-contained
# and doesn't depend on loading the full bci data
generate_layout <- function(n_trees, plot_width, plot_height, n_clusters = 10) {
  clusters <- tibble(
    cx = runif(n_clusters, 0, plot_width),
    cy = runif(n_clusters, 0, plot_height)
  )
  cluster_id <- sample(seq_len(n_clusters), n_trees, replace = TRUE)

  grid_n <- ceiling(sqrt(n_trees))
  gx <- rep(seq(0, plot_width,  length.out = grid_n), each  = grid_n)[seq_len(n_trees)]
  gy <- rep(seq(0, plot_height, length.out = grid_n), times = grid_n)[seq_len(n_trees)]

  x <- 0.55 * gx + 0.45 * clusters$cx[cluster_id] + rnorm(n_trees, 0, plot_width  * 0.015)
  y <- 0.55 * gy + 0.45 * clusters$cy[cluster_id] + rnorm(n_trees, 0, plot_height * 0.015)

  tibble(
    x = pmin(pmax(x, 0), plot_width),
    y = pmin(pmax(y, 0), plot_height)
  )
}

# ---- CORS ------------------------------------------------------------------
# Needed so index.html (opened locally / served separately) can call this API
#* @filter cors
function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  if (req$REQUEST_METHOD == "OPTIONS") {
    res$setHeader("Access-Control-Allow-Methods", "*")
    res$setHeader("Access-Control-Allow-Headers", req$HTTP_ACCESS_CONTROL_REQUEST_HEADERS)
    res$status <- 200
    return(list())
  }
  plumber::forward()
}

#* List all available species (for the species picker UI)
#* @get /species
function() {
  list(species = full_species_pool)
}

#* Add more trees (and optionally new species) to an existing forest,
#* without touching any existing trees.
#* @post /add_trees
function(n_new_trees = 100, n_new_species = 0, new_species_list = NULL,
         existing_species = NULL, next_id = 1,
         plot_width = 1000, plot_height = 500) {

  n_new_trees <- as.integer(n_new_trees)
  next_id     <- as.integer(next_id)
  plot_width  <- as.numeric(plot_width)
  plot_height <- as.numeric(plot_height)

  existing_species <- if (is.null(existing_species)) character(0) else existing_species

  # Species available for the NEW trees: whatever is already in the forest,
  # plus any newly requested species (random pick from the full ~307 pool,
  # or a specific list)
  
  new_names <- character(0)
  if (!is.null(new_species_list) && length(new_species_list) > 0) {
    new_names <- new_species_list
  } else if (as.integer(n_new_species) > 0) {
    candidates <- setdiff(full_species_pool, existing_species)
    n_pick <- min(as.integer(n_new_species), length(candidates))
    if (n_pick > 0) new_names <- sample(candidates, n_pick)
  }

  all_names <- union(existing_species, new_names)
  if (length(all_names) == 0) all_names <- sample(full_species_pool, 1)
  species_pool <- resolve_species_subset(all_names)

  layout <- generate_layout(n_new_trees, plot_width, plot_height)

  new_trees <- tibble(
    id       = seq(next_id, next_id + n_new_trees - 1),
    species  = sample(species_pool$species, n_new_trees, replace = TRUE),
    x        = layout$x,
    y        = layout$y,
    dbh_cm   = pmin(pmax(rlnorm(n_new_trees, meanlog = log(5), sdlog = 0.9), 1), 80),
    alive    = TRUE,
    age_days = 0
  )

  new_trees <- attach_species_params(new_trees, species_pool)
  new_trees <- update_allometry(new_trees)
  # Not yet shaded — will get a real value on the next /step, same as recruits.
  new_trees$shading <- 0
  new_trees$overlap <- 0

  list(
    new_trees     = new_trees,
    added_species = species_pool |> filter(!species %in% existing_species)
  )
}

#* Initialize a forest by sampling real trees from a real BCI census
#* @post /init_from_census
function(census = 1, n_trees = 2000, selected_species = NULL, plot_width = 1000, plot_height = 500) {
  census      <- as.integer(census)
  n_trees     <- as.integer(n_trees)
  plot_width  <- as.numeric(plot_width)
  plot_height <- as.numeric(plot_height)

  raw <- load_census(census)

  # Confirmed against BCI's official data documentation: A/AD/AM/AR all mean
  # "alive" (AD/AR are rare edge-case that are to be interpreated the same as A). 
  # D = dead, M = missing/ambiguous, P = not yet recruited at this census — all excluded here.
  census_data <- raw |>
    filter(status %in% c("A", "AD", "AM", "AR")) |>
    filter(!is.na(dbh) & !is.na(gx) & !is.na(gy)) |>
    mutate(dbh = dbh / 10) |>  # mm -> cm, same convention as main.R
    rename(id = treeID, x = gx, y = gy, species = sp, dbh_cm = dbh) |>
    filter(dbh_cm > 0, x >= 0, x <= plot_width, y >= 0, y <= plot_height)

  if (!is.null(selected_species) && length(selected_species) > 0) {
    census_data <- census_data |> filter(species %in% selected_species)
  }

  n_available <- nrow(census_data)
  if (n_available == 0) {
    stop(sprintf("No usable trees found in census %d after filtering.", census))
  }
  n_trees <- min(n_trees, n_available)
  
  census_data <- census_data |> mutate(.row_id = row_number())

  one_per_species <- census_data |>
    group_by(species) |>
    slice_sample(n = 1) |>
    ungroup()

  # If there are more species than the tree budget allows, we can't
  # guarantee full coverage — fall back to covering as many species as
  # n_trees allows, chosen randomly, rather than silently exceeding budget.
  if (nrow(one_per_species) > n_trees) {
    one_per_species <- one_per_species |> slice_sample(n = n_trees)
  }

  remaining_slots <- n_trees - nrow(one_per_species)
  remaining_pool <- census_data |> filter(!.row_id %in% one_per_species$.row_id)
  fill <- if (remaining_slots > 0 && nrow(remaining_pool) > 0) {
    remaining_pool |> slice_sample(n = min(remaining_slots, nrow(remaining_pool)))
  } else {
    remaining_pool[0, ]
  }

  forest <- bind_rows(one_per_species, fill) |>
    select(species, x, y, dbh_cm) |>
    mutate(alive = TRUE, age_days = 0, id = row_number())

  # attach_species_params() already falls back to global medians for any
  # species code not present in species_allometry_params.csv, so real
  # census species not in our 162 are handled 
  forest <- attach_species_params(forest, species_params)
  forest <- update_allometry(forest)
  forest <- calculate_shading(forest)

  species_in_sample <- unique(forest$species)
  species_subset <- resolve_species_subset(species_in_sample)

  list(
    forest             = forest,
    species_params     = species_subset,
    plot_width         = plot_width,
    plot_height        = plot_height,
    n_available        = n_available,
    n_sampled          = nrow(forest),
    n_species_available = length(unique(census_data$species)),
    n_species_sampled   = length(species_in_sample)
  )
}

#* Initialize a demo forest
#* @post /init
function(n_trees = 2000, n_species = 8, selected_species = NULL, plot_width = 1000, plot_height = 500) {
  n_trees     <- as.integer(n_trees)
  plot_width  <- as.numeric(plot_width)
  plot_height <- as.numeric(plot_height)

  if (!is.null(selected_species) && length(selected_species) > 0) {
    species_subset <- resolve_species_subset(selected_species)
  } else {
    n_species <- min(as.integer(n_species), length(full_species_pool))
    chosen_names <- sample(full_species_pool, n_species)
    species_subset <- resolve_species_subset(chosen_names)
  }

  layout <- generate_layout(n_trees, plot_width, plot_height)

  forest <- tibble(
    id       = seq_len(n_trees),
    species  = sample(species_subset$species, n_trees, replace = TRUE),
    x        = layout$x,
    y        = layout$y,
    dbh_cm   = pmin(pmax(rlnorm(n_trees, meanlog = log(5), sdlog = 0.9), 1), 80),
    alive    = TRUE,
    age_days = 0
  )

  forest <- attach_species_params(forest, species_subset)
  forest <- update_allometry(forest)
  forest <- calculate_shading(forest)

  list(
    forest         = forest,
    species_params = species_subset,
    plot_width     = plot_width,
    plot_height    = plot_height
  )
}

#* Advance the simulation by one real time-step
#* @post /step
function(req) {
  body <- jsonlite::fromJSON(req$postBody, simplifyDataFrame = TRUE)

  forest         <- as_tibble(body$forest)
  sim_config     <- body$sim_config
  species_subset <- as_tibble(body$species_params)

  forest <- update_allometry(forest)
  forest <- calculate_shading(forest)
  forest <- update_growth(forest, sim_config)

  mortality_result <- update_mortality(forest, sim_config)
  forest <- mortality_result$forest
  deaths <- mortality_result$deaths

  if (deaths > 0) {
    forest <- recruit_trees(forest, sim_config, species_subset)
  }
  forest <- update_allometry(forest)  # always refresh, so crown/height reflect this step's growth + any new recruits

  list(
    forest = forest,
    metrics = list(
      n_alive      = sum(forest$alive),
      mean_dbh     = mean(forest$dbh_cm[forest$alive], na.rm = TRUE),
      mean_shading = mean(forest$shading[forest$alive], na.rm = TRUE),
      deaths       = deaths
    )
  )
}
