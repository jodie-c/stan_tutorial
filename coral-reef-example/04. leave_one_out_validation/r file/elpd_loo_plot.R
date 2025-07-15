library(ggplot2)
library(dplyr)
library(bayesplot)
library(gridExtra)
library(loo)
library(patchwork) # For combining plots

data_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds"
loo_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/leave_one_out_validation"
plot_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/plots/leave_one_out_validation"

model_names <- c(
  "Complete Pooling",
  "No Pooling",
  "Partial Pooling (Centered)",
  "Partial Pooling (Non-centered)"
)

# Load loo results
model_results <- list()
for (i in 1:length(model_names)) {
  loo_file <- file.path(loo_dir, paste0("loo_", gsub(" ", "_", tolower(model_names[i])), ".rds"))
  model_results[[i]] <- readRDS(loo_file)
}

# Create a data frame with all Pareto k values
k_data <- data.frame()

for (i in 1:length(model_results)) {
  k_values <- model_results[[i]]$diagnostics$pareto_k
  
  model_k <- data.frame(
    index = 1:length(k_values),
    k = k_values,
    model = model_names[i]
  )
  
  k_data <- rbind(k_data, model_k) #index, elpd, model name(combined all the models)
}

# Create faceted plot
k_plot <- ggplot(k_data, aes(x = index, y = k, color = model)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_hline(yintercept = 0.7, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0.5, color = "orange", linetype = "dashed") +
  facet_wrap(~ model, ncol = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Pareto k Diagnostic for All Models",
    x = "Data Point Index",
    y = "Pareto k"
  ) +
  theme_minimal(base_size = 12) + #get help from chatgpt about the theme
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "lightgray", color = NA),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(plot_dir, "combined_pareto_k_plot.png"), k_plot, width = 10, height = 8)

# Combined pointwise ELPD plot for all models
elpd_data <- data.frame()

for (i in 1:length(model_results)) {
  # Get pointwise elpd values
  pointwise_elpd <- model_results[[i]]$pointwise[,"elpd_loo"]
  
  model_elpd <- data.frame(
    index = 1:length(pointwise_elpd),
    elpd = pointwise_elpd,
    model = model_names[i]
  )
  
  elpd_data <- rbind(elpd_data, model_elpd)
}

# Create faceted plot
elpd_plot <- ggplot(elpd_data, aes(x = index, y = elpd, color = model)) +
  geom_point(size = 1, alpha = 0.7) +
  facet_wrap(~ model, ncol = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Pointwise ELPD for All Models",
    x = "Data Point Index",
    y = "ELPD"
  ) +
  theme_minimal(base_size = 12) + #get help from chatgpt
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "lightgray", color = NA),
    strip.text = element_text(face = "bold")
  )

# Save the combined ELPD plot
ggsave(file.path(plot_dir, "combined_pointwise_elpd_plot.png"), elpd_plot, width = 10, height = 8)

# 3. Create a single combined diagnostic plot with both visualizations (get help from chatgpt)
combined_plot <- k_plot / elpd_plot
ggsave(file.path(plot_dir, "combined_diagnostics.png"), combined_plot, width = 12, height = 14)

cat("Combined diagnostic plots created and saved to:", plot_dir, "\n")
