# simulation_config.R

# Simulation configuration parameters

sim_config <- list(
  
  # reproducibility
  random_seed = 42,
  
  # time settings
  n_years = 40,            # total years to simulate
  n_days = 365*40,         # total days
  time_step = 365*5,
  snapshot_interval = 365*5, # save every N days
  shading_update_interval = 365*5,
  
  # # time settings
  # n_years = 40,           # 8 censuses x 5 years
  # n_days = 365*40,        # total days
  # time_step = 365*5,      # 5-year step
  # snapshot_interval = 365*5,
  # shading_update_interval = 365*5,
  
  # plot settings
  plot_width = 1000,  # meters
  plot_height = 500,  # meters
  
  # spatial settings
  quad_size = 46,     # quadrant size for shading
  
  # initial conditions
  inital_census = "census8", # which census to start
  
  # growth parameters
  growth = list(
    beta_0 = -0.935,
    beta_D = 0.369,
    beta_S = -0.453
  ),
  
  # mortality parameters
  mortality = list(
    beta_0 = -1.319,
    beta_D = -0.303,
    beta_S = -.172
  ),
  
  # recruitment parameters
  recruitment = list(
    fecundity_power = 2, # DBH exp controlling fec
    sigma_seed = 30, # seed dispersal spread (10)
    neighbor_radius = 100, # max distance to search for parents
    recruit_dbh = 1 # initial DBH of recruit is 1cm
  ),
  
  
  # output settings
  save_snapshots = TRUE,
  save_metrics = TRUE,
  verbose = TRUE
  
)

