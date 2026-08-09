# allometry.R

# calculate tree height and crown area from DBH

attach_species_params <- function(forest, species_params){
  
  
  # compute global medians for missing species
  global_crown_a <- median(species_params$a, na.rm = TRUE)
  global_crown_b <- median(species_params$b, na.rm = TRUE)
  
  global_height_a <- median(species_params$ta, na.rm = TRUE)
  global_height_b <- median(species_params$tb, na.rm = TRUE)
  global_height_k <- median(species_params$tk, na.rm = TRUE)
  
  
  forest <- forest |>
    left_join(
      species_params |>
        select(
          species, 
          crown_a = a, 
          crown_b = b,
          height_a = ta, 
          height_b = tb, 
          height_k = tk
        ),
      by = "species"
    ) |>
    mutate(
      crown_a = ifelse(is.na(crown_a), global_crown_a, crown_a),
      crown_b = ifelse(is.na(crown_b), global_crown_b, crown_b),
      
      height_a = ifelse(is.na(height_a), global_height_a, height_a),
      height_b = ifelse(is.na(height_b), global_height_b, height_b),
      height_k = ifelse(is.na(height_k), global_height_k, height_k)
    )
}

update_allometry <- function(forest){
  
  # calculate height, crown area, and radius (size): 
  # H = a * D^b / (k + D^b)
  # A = a * D^b
  # size = sqrt(A / pi)
  
  forest <- forest |>
    mutate(
      # height
      height_m = (height_a*dbh_cm^height_b) / (height_k+dbh_cm^height_b),
      
      # crown area
      crown_area_m2 = crown_a*dbh_cm^crown_b,
      
      # crown radius aka 'size' (for overlap calc)
      crown_radius_m = sqrt(crown_area_m2/pi)
    )
  
  return(forest)
}


