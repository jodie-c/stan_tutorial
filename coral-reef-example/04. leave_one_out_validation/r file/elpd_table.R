library(rstan)
library(loo)
library(dplyr)

data_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds"
output_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/leave_one_out_validation"

fit_files <- c(
  "stan_fit_complete_pooling.rds",
  "stan_fit_no_pooling.rds", 
  "stan_fit_partial_centered.rds",
  "stan_fit_partial_noncentered.rds"
)

model_names <- c(
  "Complete Pooling",
  "No Pooling",
  "Partial Pooling (Centered)",
  "Partial Pooling (Non-centered)"
)

# Get the data
data_file <- file.path(data_dir, "stan_data_complete_pooling.rds")
model_data <- readRDS(data_file)

n_sites <- model_data$N_sites
n_years <- model_data$N_years
n_data_points <- sum(model_data$Y_obs_present)
t_grid <- model_data$t_gd
time_points <- model_data$ts

model_results <- list()
log_lik_data <- list()

# Process each model
for (i in 1:length(fit_files)) {
  
  fit <- readRDS(file.path(data_dir, fit_files[i]))
  
  # Extract samples
  y_fit_tot <- rstan::extract(fit, pars = "y_fit_tot")$y_fit_tot
  sigma_samples <- rstan::extract(fit, pars = "sigma")$sigma
  
  n_samples <- dim(y_fit_tot)[1]
  
  # Calculate log-likelihood matrix of data point under sample s(full posterior distribution)
  log_lik <- matrix(0, nrow = n_samples, ncol = n_data_points) # rows: posterior samples, columns: data points, value: log-likelihood
  # For each posterior sample s
  for (s in 1:n_samples) {
    point_idx <- 1
    
    # For each site r
    for (r in 1:n_sites) {
      # For each time point t
      for (t in 1:n_years) {
        if (model_data$Y_obs_present[t, r] == 1) { #skip if the site is not present in the year(for the sparse data)
          # match prediction time to the closest time in the dense grid
          closest_t_idx <- which.min(abs(t_grid - time_points[t]))
          
          # Get predicted and observed coral cover(actual parameter values)
          y_pred <- y_fit_tot[s, r, closest_t_idx] #sample s, site r, the matched time point
          y_obs <- model_data$Y_obs[t, r] #observed coral cover
          
          # Handle different sigma structures across models, complete pooling one sigma for all sites, other pooling methods have different sigmas for differnt sites
          if (length(dim(sigma_samples)) > 1) {
            # Site-specific sigma (no pooling, partial pooling models)
            sigma_r <- sigma_samples[s, r] 
          } else {
            # Shared sigma (complete pooling model)
            sigma_r <- sigma_samples[s]
          }
          
          # likelihood of data point d under the parameter values from posterior sample s
          log_lik[s, point_idx] <- dnorm(y_obs, mean = y_pred, sd = sigma_r, log = TRUE) 
          
          point_idx <- point_idx + 1 
        }
      }
    }
  }
  
  # Save log-likelihood matrix
  log_lik_data[[i]] <- log_lik
  
  # Compute LOO
  model_results[[i]] <- loo(log_lik)
  print(model_results[[i]])
  
  # Save individual LOO results
  saveRDS(model_results[[i]], file.path(output_dir, paste0("loo_", gsub(" ", "_", tolower(model_names[i])), ".rds")))
}

# Save all log-likelihood matrices
saveRDS(log_lik_data, file.path(output_dir, "all_log_lik.rds"))

# Compare models
if (length(model_results) > 1) {
  cat("\nModel comparison:\n")
  model_comparison <- loo_compare(model_results)
  print(model_comparison)
  
  # Save comparison results
  saveRDS(model_comparison, file.path(output_dir, "loo_comparison.rds"))
  
  # Create a summary table
  elpd_values <- numeric(length(model_results))
  se_values <- numeric(length(model_results))
  ploo_values <- numeric(length(model_results))
  looic_values <- numeric(length(model_results))
  
  # Extract values from each model's results
  for (i in 1:length(model_results)) {
    elpd_values[i] <- model_results[[i]]$elpd_loo    # ELPD value
    se_values[i] <- model_results[[i]]$se_elpd_loo   # ELPD standard error
    ploo_values[i] <- model_results[[i]]$p_loo       # P_LOO value
    looic_values[i] <- model_results[[i]]$looic      # LOOIC value
  }      
  
  model_summary <- data.frame(
    Model = model_names,
    ELPD = elpd_values,
    SE = se_values,
    P_LOO = ploo_values,
    LOOIC = looic_values
  )
  
  # Sort by ELPD value (higher is better)
  model_summary <- model_summary %>% arrange(desc(ELPD))
  print(model_summary)
  
  write.csv(model_summary, file.path(output_dir, "loo_summary.csv"), row.names = FALSE)
  
 
 
}

