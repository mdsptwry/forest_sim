# shading.R

# calculate shading for all trees

calculate_shading <- function(forest, quad_size = 46){
  
  # remove old shading columns
  # the reason we do this is because the initial forest doesn't have shading
  # but the subsequent ones will, so we remove the old ones
  # and replace them with the new shading
  forest <- forest |>
    select(
      -any_of(c("shading", "overlap"))
    )
  
  # forest prep 
  # only compute shading for living trees
  # Create size variable which is crown_radius_m
  forest_prep <- forest |>
    filter(alive) |>
    mutate(size = crown_radius_m) |>
    select(id, x, y, size, height_m)
  
  
  overlap_results <- calculate_overlaps_knn(forest_prep)
  
  # extract shading
  forest <- forest |>
    left_join(
      overlap_results |>
        select(id, overlap, shading),
      by = "id"
    ) |>
    mutate(
      shading = ifelse(is.na(shading), 0, shading),
      overlap = ifelse(is.na(overlap), 0, overlap)
    )
  
  return(forest)
  
}

