library(rstan)
library(dplyr)
library(tidyr)

data_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds"

#step 1: load the models
models <- list(
  list(
    name = "Complete Pooling",
    fit_file = file.path(data_dir, "stan_fit_complete_pooling_normal.rds"),
    data_file = file.path(data_dir, "stan_data_complete_pooling_normal.rds")
  ),
  list(
    name = "No Pooling", 
    fit_file = file.path(data_dir, "stan_fit_no_pooling.rds"),
    data_file = file.path(data_dir, "stan_data_no_pooling.rds")
  ),
  list(
    name = "Partial Pooling Centered",
    fit_file = file.path(data_dir, "stan_fit_partial_centered.rds"),
    data_file = file.path(data_dir, "stan_data_partial_centered.rds")
  ),
  list(
    name = "Partial Pooling Non-Centered",
    fit_file = file.path(data_dir, "stan_fit_partial_noncentered.rds"),
    data_file = file.path(data_dir, "stan_data_partial_noncentered.rds")
  )
)

low_threshold <- 0.05
high_threshold <- 0.95

#step 2: calculate the 7 p-values
model_results <- list()

for (i in 1:length(models)) {
  model_info <- models[[i]]
  fit <- readRDS(model_info$fit_file)
  stan_data <- readRDS(model_info$data_file)
  
  obs_matrix <- stan_data$Y_obs # [10, 45] time x site matrix of observed coral cover
  obs_present <- stan_data$Y_obs_present # [10, 45] binary matrix indicating which observations exist
  n_sites <- stan_data$N_sites #number of reef sites (45)
  n_years <- stan_data$N_years #number of time points (10)
  
  rep_array <- rstan::extract(fit, pars = "Y_rep")$Y_rep # [4000, 10, 45] iterations x time x site posterior predictive samples
  obs_vector <- as.vector(obs_matrix) # [450] flattened vector of all observed values (complete data)
  n_observed <- length(obs_vector) # total number of data points (450)
  n_iter <- dim(rep_array)[1] # number of iterations (4000)
  rep_matrix <- matrix(rep_array, nrow = n_iter, ncol = n_years * n_sites) # [4000, 450] all posterior samples reshaped
  
  rep_matrix[rep_matrix == -999] <- NA # replace missing value indicators with NA
  
    #mean
    # [450]observed data and then calculate the mean value, one single number
    obs_mean <- mean(obs_vector, na.rm = TRUE) 
    # each row is one posterior sample, and then calculate the mean value, [4000] vector
    rep_mean <- apply(rep_matrix, 1, mean, na.rm = TRUE) 
    #how many out of 4000 posterior samples are greater than the observed mean?
    p_mean <- mean(rep_mean >= obs_mean)

    obs_sd <- sd(obs_vector, na.rm = TRUE)
    rep_sd <- apply(rep_matrix, 1, sd, na.rm = TRUE)
    p_sd <- mean(rep_sd >= obs_sd)
    
    obs_min <- min(obs_vector, na.rm = TRUE)
    rep_min <- apply(rep_matrix, 1, min, na.rm = TRUE)
    p_min <- mean(rep_min >= obs_min)
    
    obs_max <- max(obs_vector, na.rm = TRUE)
    rep_max <- apply(rep_matrix, 1, max, na.rm = TRUE)
    p_max <- mean(rep_max >= obs_max)
    
    obs_q25 <- quantile(obs_vector, 0.25, na.rm = TRUE)
    rep_q25 <- apply(rep_matrix, 1, quantile, 0.25, na.rm = TRUE)
    p_q25 <- mean(rep_q25 >= obs_q25)
    
    obs_q75 <- quantile(obs_vector, 0.75, na.rm = TRUE)
    rep_q75 <- apply(rep_matrix, 1, quantile, 0.75, na.rm = TRUE)
    p_q75 <- mean(rep_q75 >= obs_q75)
    
    obs_median <- median(obs_vector, na.rm = TRUE)
    rep_median <- apply(rep_matrix, 1, median, na.rm = TRUE)
    p_median <- mean(rep_median >= obs_median)
  
  p_values <- c(p_mean, p_sd, p_min, p_max, p_q25, p_q75, p_median) # [7] vector of all calculated p-values
  valid_p_values <- na.omit(p_values) # [7] vector with any NA values removed
  
        #percentage of extreme p-values
   pct_extreme <- mean(valid_p_values < low_threshold | valid_p_values > high_threshold) * 100 
  
   # Calculate average distance from 0.5 for all p-values
   avg_distance <- mean(abs(valid_p_values - 0.5), na.rm = TRUE)
  
  model_results[[i]] <- data.frame(
    Model = model_info$name,
    p_mean,
    p_sd,
    p_min,
    p_max,
    p_q25,
    p_q75,
    p_median,
    pct_extreme,
    avg_distance,
    check.names = FALSE
  )
}

# step3: display the results
results_df <- bind_rows(model_results) # data frame: [4 x 17] combined results from all models
print.data.frame(as.data.frame(results_df), digits = 3, na.print = "NA", row.names = FALSE)


cat("SUMMARY TABLE: Average Distance from 0.5 and Extreme P-values\n")

final_table <- results_df %>% # data frame: [4 x 4] 
  rowwise() %>%
  mutate(
    avg_pvalue = mean(c(p_mean, p_sd, p_min, p_max, 
                       p_q25, p_q75, p_median), na.rm = TRUE)
  ) %>%
  select(Model, avg_distance, pct_extreme, avg_pvalue) %>%
  arrange(avg_distance) %>%
  mutate(
    avg_distance = round(avg_distance, 3),
    pct_extreme = round(pct_extreme, 3),
    avg_pvalue = round(avg_pvalue, 3)
  )

cat("\nFINAL TABLE:\n")
print(final_table, row.names = FALSE)
