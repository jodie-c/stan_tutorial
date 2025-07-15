library(rstan)

model_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds"
complete_pooling <- readRDS(file.path(model_dir, "stan_fit_complete_pooling.rds"))
no_pooling <- readRDS(file.path(model_dir, "stan_fit_no_pooling.rds"))
partial_centered <- readRDS(file.path(model_dir, "stan_fit_partial_centered.rds"))
partial_noncentered <- readRDS(file.path(model_dir, "stan_fit_partial_noncentered.rds"))

# Complete Pooling
model_summary <- summary(complete_pooling)
all_params <- rownames(model_summary$summary)
core_params <- all_params[!grepl("y_fit_tot|Y_rep|log_lik|lp__", all_params)]
param_count <- length(core_params)
divergent_count <- sum(get_divergent_iterations(complete_pooling))
rhat_mean <- mean(model_summary$summary[core_params, "Rhat"], na.rm = TRUE)

# No Pooling
model_summary <- summary(no_pooling)
all_params <- rownames(model_summary$summary)
core_params <- all_params[!grepl("y_fit_tot|Y_rep|log_lik|lp__", all_params)]
param_count_np <- length(core_params)
divergent_count_np <- sum(get_divergent_iterations(no_pooling))
rhat_mean_np <- mean(model_summary$summary[core_params, "Rhat"], na.rm = TRUE)

# Partial Pooling Centered
model_summary <- summary(partial_centered)
all_params <- rownames(model_summary$summary)
core_params <- all_params[!grepl("y_fit_tot|Y_rep|log_lik|lp__", all_params)]
param_count_pc <- length(core_params)
divergent_count_pc <- sum(get_divergent_iterations(partial_centered))
rhat_mean_pc <- mean(model_summary$summary[core_params, "Rhat"], na.rm = TRUE)

# Partial Pooling Non-Centered
model_summary <- summary(partial_noncentered)
all_params <- rownames(model_summary$summary)
core_params <- all_params[!grepl("y_fit_tot|Y_rep|log_lik|lp__", all_params)]
param_count_pnc <- length(core_params)
divergent_count_pnc <- sum(get_divergent_iterations(partial_noncentered))
rhat_mean_pnc <- mean(model_summary$summary[core_params, "Rhat"], na.rm = TRUE)

# Print results
cat("Complete Pooling:", param_count, "parameters,", divergent_count, "divergent,", round(rhat_mean, 4), "rhat\n")
cat("No Pooling:", param_count_np, "parameters,", divergent_count_np, "divergent,", round(rhat_mean_np, 4), "rhat\n")
cat("Partial Pooling Centered:", param_count_pc, "parameters,", divergent_count_pc, "divergent,", round(rhat_mean_pc, 4), "rhat\n")
cat("Partial Pooling Non-Centered:", param_count_pnc, "parameters,", divergent_count_pnc, "divergent,", round(rhat_mean_pnc, 4), "rhat\n")
