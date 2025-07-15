library(rstan)
library(shinystan)

data_dir <- '/Users/chengu/Documents/master thesis/Bayesian-inference-of-coral-bleaching-dynamics/01. pooling_methods/rds'

launch_model_viewer <- function(model_name) {
  model_file <- file.path(data_dir, paste0("stan_fit_", model_name, ".rds"))
  
  # Load the model
  model <- readRDS(model_file)
  
  # Convert to shinystan object
  shiny_object <- as.shinystan(model)
  
  launch_shinystan(shiny_object, launch.browser = FALSE)
}

launch_model_viewer("partial_noncentered")