# mortality.R

# determine which trees die

update_mortality <- function(forest, sim_config){
  
  # get params
  beta_0 <- sim_config$mortality$beta_0
  beta_D <- sim_config$mortality$beta_D
  beta_S <- sim_config$mortality$beta_S
  
  # calculate daily death probability
  forest <- forest |>
    mutate(
      # log-odds of death
      log_odds_death = ifelse(
        alive & dbh_cm > 0,
        beta_0 + beta_D*log(dbh_cm) + beta_S*shading,
        0
      ),
      
      # iteration, change formula to +
      
      # probability of death (5-year)
      p_death_5yr = plogis(log_odds_death),
      
      # convert to daily probability
      # if 5-year survival = (1-p)^(5*365)
      # then daily p = 1-(1-p_5yr)^(1/(5*365))
      
      #p_death_daily = 1-(1-p_death_5yr)^(1/(5*365)),
      p_death_step = 1-(1-p_death_5yr)^(sim_config$time_step/(5*365)),
      
      # stochastic death
      dies = alive & (runif(n()) < p_death_step),
    )
  
  deaths_this_step <- sum(forest$dies, na.rm = TRUE)
  
  forest <- forest |>
    mutate(
      alive = alive & !dies
    ) |>
    select(
      -log_odds_death, 
      -p_death_5yr, 
      -p_death_step)
  
  return(list(
    forest = forest,
    deaths = deaths_this_step
  ))
}

