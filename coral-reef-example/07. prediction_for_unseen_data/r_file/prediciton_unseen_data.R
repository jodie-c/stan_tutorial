FORGETTING_FACTOR <- 0.9

library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(deSolve)
library(purrr)

data_dir <- '/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds'
raw_data_file <- '/Users/chengu/Documents/master thesis/new_reef_data.csv'
fit_file <- file.path(data_dir, "stan_fit_no_pooling.rds")
data_file <- file.path(data_dir, "stan_data_no_pooling.rds")
plot_dir <- file.path('/Users/chengu/Documents/master thesis/final_final_final_code/plots/prediction_for_unseen_data', "validation_plots_from_fit_model")

site_colors <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A")
START_YEAR <- 1995
TRAIN_END <- 2004
VALID_END <- 2007

fit <- readRDS(fit_file)
stan_data <- readRDS(data_file)

reef_names <- readRDS(file.path(data_dir, "reef_id_mapping.rds"))$REEF_NAME

# Create site mapping
site_map <- data.frame(
  site_id = 1:45,
  reef_name = rep(reef_names, each = 3),
  site_number = rep(1:3, 15),
  reef_id = rep(1:15, each = 3)
)

# Load and process raw data
raw_data <- read.csv(raw_data_file, sep = ';', header = TRUE, stringsAsFactors = FALSE) %>%
  select(REEF_NAME, SITE_NO, SAMPLE_DATE, GROUP_CODE, COVER)

# Process observed data
reef_data <- raw_data %>%
  filter(GROUP_CODE == 'Hard Coral') %>%
  mutate(
    YEAR = as.numeric(substr(SAMPLE_DATE, 1, 4)),
    Cover = as.numeric(gsub(',', '.', COVER)) / 100
  ) %>%
  filter(YEAR >= START_YEAR, YEAR <= VALID_END) %>%
  group_by(REEF_NAME, YEAR, SITE_NO) %>%
  summarize(Cover = mean(Cover, na.rm = TRUE), .groups = 'drop') %>%
  inner_join(site_map, by = c("REEF_NAME" = "reef_name", "SITE_NO" = "site_number")) %>%
  mutate(Period = ifelse(YEAR <= TRAIN_END, "Training", "Validation"))

# Extract parameters
# params$theta shape: [n_samples, 4, n_sites], a 3D array: n_samples x 4 x n_sites
# params$y0_init_C shape: [n_samples, n_sites]
params <- rstan::extract(fit, pars = c("theta", "y0_init_C"))

# Generate predictions
prediction_times <- seq(0, VALID_END - START_YEAR, by = 0.1)
set.seed(123)
n_samples <- min(1000, dim(params$theta)[1])
sample_indices <- sample(n_samples)

# Solve ODEs
# all_trajectories: sample, site_id, time, total_cover
all_trajectories <- list()

for (i in seq_along(sample_indices)) {
  sample_idx <- sample_indices[i]
  
  for (s in 1:stan_data$N_sites) {
    # theta_s shape: [4] - vector with alpha, beta, gamma, mu parameters
    theta_s <- pmax(1e-6, params$theta[sample_idx, , s])
    N0 <- params$y0_init_C[sample_idx, s]
    # y0 shape: [2] - vector with initial conditions for C and B 
    y0 <- c(C = pmax(1e-6, pmin(1-1e-6, N0^2)), 
            B = pmax(1e-6, pmin(1-1e-6, N0 * (1-N0))))
    
    # ODE model
    # parms[1] = alpha, parms[2] = beta 
    # parms[3] = gamma, parms[4] = mu 
    coral_ode <- function(t, y, parms) {
      C <- pmax(0, pmin(1, y[1]))
      B <- pmax(0, pmin(1, y[2]))
      available <- pmax(0, 1 - C - B)
      dCdt <- parms[1] * C * available - parms[2] * C + parms[3] * B
      dBdt <- parms[2] * C - parms[3] * B - parms[4] * B
      list(c(dCdt, dBdt))
    }
    
    # Solve the ODE system
    # ode_result: time points, C, B
    ode_result <- deSolve::ode(y = y0, times = prediction_times, func = coral_ode, parms = theta_s)
    
    all_trajectories[[length(all_trajectories) + 1]] <- data.frame(
      sample = sample_idx, 
      site_id = s, 
      time = ode_result[,1],
      total_cover = pmin(1, pmax(0, rowSums(ode_result[,2:3])))
    )
  }
}

prediction_summary <- bind_rows(all_trajectories) %>%
  group_by(site_id, time) %>%
  summarize(
    lower_ci = quantile(total_cover, 0.025, na.rm = TRUE),
    median_cover = quantile(total_cover, 0.5, na.rm = TRUE),
    upper_ci = quantile(total_cover, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(YEAR = time + START_YEAR) %>% 
  left_join(site_map, by = c("site_id" = "site_id"))

# Make sure reef names appear in the same order in both datasets (get help from chapgpt these two lines)

# prediction_summary$reef_name: A character vector with reef names for each row in the dataframe
# levels = reef_names: Sets the order of categories using our predefined reef_names vector
prediction_summary$reef_name <- factor(prediction_summary$reef_name, levels = reef_names)
reef_data$REEF_NAME <- factor(reef_data$REEF_NAME, levels = reef_names)

# Create validation plot (get help from chatgpt about the theme of the plot)
validation_plot <- ggplot() +
  geom_ribbon(data = prediction_summary, 
              aes(x = YEAR, ymin = lower_ci, ymax = upper_ci, 
                  group = factor(site_number), fill = factor(site_number)), 
              alpha = 0.3) +
  geom_line(data = prediction_summary, 
            aes(x = YEAR, y = median_cover, 
                group = factor(site_number), color = factor(site_number)), 
            linewidth = 1) +
  geom_point(data = reef_data, 
             aes(x = YEAR, y = Cover, color = factor(SITE_NO), shape = Period), 
             size = 2.5, alpha = 0.8) +
  geom_vline(xintercept = TRAIN_END + 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  scale_color_manual(values = site_colors, name = "Site") +
  scale_fill_manual(values = site_colors, name = "Site") +
  scale_shape_manual(name = "Period", values = c("Training" = 16, "Validation" = 17)) +
  facet_wrap(~ reef_name, scales = "free_y", ncol = 3) +
  labs(title = "Predicted vs. Observed Coral Cover (1995-2007)",
       subtitle = "Lines are median predictions, ribbons are 95% credible intervals. Circles are training data, triangles are validation data.",
       x = "Year", y = "Total Coral Cover") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

output_file <- file.path(plot_dir, "validation_plot_no_pooling.pdf")
ggsave(output_file, plot = validation_plot, width = 15, height = 20, device = "pdf")
