# intialization.R

# initialize forest from census data

initialize_forest <- function(census_data, species_params, sim_config){
  
  # start with census data
  forest <- census_data |>
    select(id, species, x, y, dbh_cm) |>
    mutate(
      alive = TRUE,
      age_days = 0
    )
  
  # attach species params
  forest <- attach_species_params(forest, species_params)
  
  # update allometry
  forest <- update_allometry(forest)
  
  # calculate shading
  forest <- calculate_shading(forest, quad_size = sim_config$quad_size)
  
  return(forest)
  
}

