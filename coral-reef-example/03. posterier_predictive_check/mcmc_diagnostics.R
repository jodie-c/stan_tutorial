library(rstan)
library(dplyr)

data_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds"

models <- list(
  list(
    name = "Complete Pooling",
    fit_file = file.path(data_dir, "stan_fit_complete_pooling.rds")
  ),
  list(
    name = "No Pooling", 
    fit_file = file.path(data_dir, "stan_fit_no_pooling.rds")
  ),
  list(
    name = "Partial Pooling Centered",
    fit_file = file.path(data_dir, "stan_fit_partial_centered.rds")
  ),
  list(
    name = "Partial Pooling Non-Centered",
    fit_file = file.path(data_dir, "stan_fit_partial_noncentered.rds")
  )
)


# Create results table
model_stats <- data.frame(
  Model = character(),
  Parameters = integer(),
  Divergent_Transitions = integer(),
  Average_Rhat = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:length(models)) {
  model_info <- models[[i]]
  
  fit <- readRDS(model_info$fit_file)
  
  # Get basic info
  summary_fit <- summary(fit)
  param_names <- rownames(summary_fit$summary)
  
  # Count parameters (exclude log posterior)
  total_params <- sum(!grepl("lp__", param_names))
  
  # Get divergent transitions
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  divergent_count <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
  
  # Get average Rhat (exclude log posterior)
  rhat_values <- summary_fit$summary[, "Rhat"]
  # Filter out log posterior and other diagnostics
  clean_rhat <- rhat_values[!grepl("lp__", param_names)]
  clean_rhat <- clean_rhat[!is.na(clean_rhat) & is.finite(clean_rhat)]
  avg_rhat <- mean(clean_rhat)
  
  # Add to results
  model_stats <- rbind(model_stats, data.frame(
    Model = model_info$name,
    Parameters = total_params,
    Divergent_Transitions = divergent_count,
    Average_Rhat = round(avg_rhat, 4),
    stringsAsFactors = FALSE
  ))
}

# Display results(get help from chatgpt)
cat("MCMC DIAGNOSTICS SUMMARY:\n")
cat("=========================\n\n")

cat(sprintf("%-30s %12s %20s %15s\n", 
            "Model", "Parameters", "Divergent Transitions", "Average Rhat"))
cat(sprintf("%-30s %12s %20s %15s\n", 
            "-----", "----------", "-------------------", "------------"))

for(i in 1:nrow(model_stats)) {
  row <- model_stats[i,]
  cat(sprintf("%-30s %12d %20d %15.4f\n", 
              row$Model, row$Parameters, row$Divergent_Transitions, row$Average_Rhat))
}

cat("\nDiagnostic Explanations:\n")
cat("Parameters:")
cat("Divergent Transitions:")
cat("Average Rhat:")

# Save results
write.csv(model_stats, "mcmc_diagnostics.csv", row.names = FALSE)
cat("\nResults saved to: mcmc_diagnostics.csv\n")
