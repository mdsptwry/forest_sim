# main.R

# forest simulation - main script

cat("Forest Dynamics Simulation\n")

# load libraries
library(tidyverse)
library(furrr)
library(ggforce)
library(patchwork)
library(dplyr)
library(ggplot2)


# set up parallel processing
plan(multisession(workers = availableCores()-2))


# load configuration
source("config/simulation_config.R")

# set random seed for reproducibility
set.seed(sim_config$random_seed)
options(future.seed = sim_config$random_seed)

# load all functions/files
source("R/allometry.R")
source("R/spatial_functions.R")
source("R/shading.R")
source("R/initialization.R")
source("R/growth.R")
source("R/mortality.R")
source("R/recruitment.R")
source("R/utilities.R")

cat("Load species parameters...\n")

species_params <- read.csv("data/species_allometry_params.csv")

# initialize forest

cat("Initializing forest...\n")
load("data/bci.tree2.rdata")

census_data <- bci.tree2 |>
  filter(status %in% c("A", "AD", "AM", "AR")) |>   # alive-only (confirmed BCI status codes)
  filter(!is.na(dbh) & !is.na(gx) & !is.na(gy)) |>
  mutate(
    dbh = dbh/10, # convert mm to cm
  ) |> rename(
    id = treeID,
    x = gx, 
    y = gy, 
    species = sp,
    dbh_cm = dbh
  )

inital_census <- census_data

forest <- initialize_forest(inital_census, species_params, sim_config)

cat ("Initial forest state:\n")
cat("Trees: ", sum(forest$alive), "\n")
cat("Species: ", length(unique(forest$species)), "\n\n")

save_snapshot(forest, day=0)


# Main simulation loop
cat("Starting Simulation\n")
cat("Duration: ", sim_config$n_years, "years (", sim_config$n_days, "days)\n")


# initialize metrics storage
metrics_list <- list()
metrics_list[[1]] <- calculate_metrics(forest, day=0)

# start time
start_time <- Sys.time()
cat("Start time: ", start_time, "\n")

# interval updates
next_snapshot <- sim_config$snapshot_interval
next_shading_update <- sim_config$shading_update_interval

# main loop
for (day in seq(sim_config$time_step, sim_config$n_days, by = sim_config$time_step)){
  
  # 1. Update allometry
  forest <- update_allometry(forest)
  cat("Allometry Updated!\n")
  
  # 2. Calculate shading
  if (day >= next_shading_update || day == 1){
    forest <- calculate_shading(forest, quad_size = sim_config$quad_size)
    cat("Shading Calculated!\n")
    next_shading_update <- next_shading_update+sim_config$shading_update_interval
  }
  
  # 3. Growth
  forest <- update_growth(forest, sim_config)
  cat("Growth Updated!\n")
  
  # 4. Mortality
  mortality_result <- update_mortality(forest, sim_config)
  forest <- mortality_result$forest
  deaths_this_step <- mortality_result$deaths
  cat("Mortality Updated!\n")
  cat("Deaths this step: ", deaths_this_step, "\n")
  
  # 5. Recruitment (only if deaths occurred)
  n_dead <- sum(!forest$alive)
  
  if (n_dead > 0){
    forest <- recruit_trees(forest, sim_config, species_params)
    forest <- update_allometry(forest)
    cat("New Recruits Added!\n")
  }
  
  
  # 6. Save snapshots
  if (day >= next_snapshot){
    
    if (sim_config$save_snapshots)
      save_snapshot(forest, day)
    
    metrics_list[[length(metrics_list)+1]] <-
      calculate_metrics(forest, day)
    
    if(sim_config$verbose)
      print_progress(forest, day, sim_config$n_days)
    
    next_snapshot <- next_snapshot + sim_config$snapshot_interval
  }
}


# Save Results

saveRDS(forest, "output/final_forest_state.rds")

metrics_df <- bind_rows(metrics_list)

write.csv(
  metrics_df,
  "output/metrics/simulation_metrics.csv",
  row.names = FALSE
)

# plots

p1 <- ggplot(metrics_df, aes(day/365, n_alive))+
  geom_line()+
  labs(x="Years", y="Living Trees")

p2 <- ggplot(metrics_df, aes(day/365, mean_dbh))+
  geom_line()+
  labs(x="Years", y="Mean DBH")

ggsave("output/plots/population.png", p1, width = 8, height = 5)
ggsave("output/plots/dbh.png", p2, width = 8, height = 5)

final_plot <- plot_forest(forest)
ggsave("output/plots/final_forest.png", final_plot, width = 12, height = 6)


end_time <- Sys.time()
cat("End Time: ", end_time, "\n")

elapsed <- difftime(end_time, start_time, units = "mins")

cat("\n Simulation Finished in", round(elapsed, 2), "minutes\n")

