# utilities.R

# helper functions for simulation output, monitoring, and visualization

# save forest snapshot

save_snapshot <- function(forest, day, output_dir = "output/snapshots"){
  
  filename <- file.path(output_dir, paste0("forest_day_", day, ".rds"))
  
  saveRDS(forest, filename)
}


# calculate forest metrics

calculate_metrics <- function(forest, day){
  
  tibble::tibble(
    day = day,
    n_alive = sum(forest$alive),
    n_dead = sum(!forest$alive),
    mean_dbh = mean(forest$dbh_cm[forest$alive], na.rm = TRUE),
    mean_shading = mean(forest$shading[forest$alive], na.rm = TRUE),
    n_species = length(unique(forest$species[forest$alive]))
  )
}


# print simulation progress

print_progress <- function(forest, day, total_days){
  
  pct_complete <- round(100*day/total_days, 1)
  
  cat(sprintf(
    "Day %d/%d (%s%%) | Alive: %d | Mean DBH: %.2f cm | Mean Shading: %.3f\n",
    day,
    total_days,
    pct_complete,
    sum(forest$alive),
    mean(forest$dbh_cm[forest$alive], na.rm = TRUE),
    mean(forest$shading[forest$alive], na.rm = TRUE)
  ))
}


# plot forest state
plot_forest <- function(forest){
  library(ggforce)
  
  forest_alive <- forest |> filter(alive)
  
  ggplot(forest_alive, aes(x = x, y = y, color = species))+
    geom_point(size = 0.3, alpha = 0.5)+
    coord_fixed(xlim = c(0, 1000), ylim = c(0, 500))+
    labs(title = "Forest", x= "X (m)", y="Y (m)", color = "Species")+
    theme_minimal()+
    guides(color = "none")
}


