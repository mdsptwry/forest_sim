# growth.R

# update DBH based on growth model

update_growth <- function(forest, sim_config){
  
  # get params
  beta_0 <- sim_config$growth$beta_0
  beta_D <- sim_config$growth$beta_D
  beta_S <- sim_config$growth$beta_S
  
  # calculate growth for live trees
  forest <- forest |>
    mutate(
      log_delta_D_5yr = ifelse(
        alive,
        beta_0 + beta_D*log(dbh_cm) + beta_S*shading, # removed log
        0
      ),
      delta_D_5yr = exp(log_delta_D_5yr),
      
      # convert to daily growth
      #daily_growth = delta_D_5yr / (5*365),
      
      growth_step = delta_D_5yr * (sim_config$time_step /(5*365)),
      
      # update DBH
      dbh_cm = ifelse(
        alive, 
        dbh_cm+growth_step, # not daily_growth, 
        dbh_cm),
      
      # increment age (is it really needed? idk)
      
      age_days = ifelse(alive, age_days+1, age_days)
    ) |>
    select(-log_delta_D_5yr, -delta_D_5yr, -growth_step)
  
  return(forest)
  
}

