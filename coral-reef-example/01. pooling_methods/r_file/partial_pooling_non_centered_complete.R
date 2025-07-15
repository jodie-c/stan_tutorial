library(rstan)
library(dplyr)
library(tidyr)
library(magrittr)

data_file <- "/Users/chengu/Documents/master thesis/new_reef_data.csv"
model_file <- "/Users/chengu/Documents/master thesis/final_final_final_code/stan_file_for_pooling_method/partial_pooling_non_centered_complete_0703.stan"
output_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data_output_file <- file.path(output_dir, "stan_data_partial_noncentered.rds")
fit_output_file <- file.path(output_dir, "stan_fit_partial_noncentered.rds")

options(mc.cores = parallel::detectCores()) 
rstan_options(auto_write = TRUE)             

# Step 1: Load and clean data
data_raw_full <- read.table(
  data_file, 
  header = TRUE, 
  sep = ";", 
  stringsAsFactors = FALSE, 
  quote = ""
) 

# Select relevant columns and process data types
data_raw <- data_raw_full %>%
  select(
    REEF_NAME, SAMPLE_DATE, SITE_NO, 
    GROUP_CODE, COVER, LATITUDE, LONGITUDE
  ) %>% 
  mutate(
    COVER = as.numeric(gsub("[^0-9.]", "", COVER)),
    SAMPLE_DATE = as.Date(SAMPLE_DATE, format = "%Y-%m-%d"),
    YEAR = as.numeric(format(SAMPLE_DATE, "%Y"))
  )

# Calculate coral cover per site
site_coral <- data_raw %>%
  # Filter for Hard Coral data only
  filter(GROUP_CODE %in% c("Hard Coral")) %>% 
  # Group by reef, year, and site
  group_by(REEF_NAME, YEAR, SITE_NO) %>%
  # Calculate total coral cover (divide by 100 to get proportion)
  summarize(site_coral = sum(COVER) / 100, .groups = "drop")

# Step 2: Select 15 reefs with complete data 1995-2004
site_coral_filtered <- site_coral %>% 
  filter(YEAR >= 1995, YEAR <= 2004)

# Count number of sites per reef per year
sites_per_year <- site_coral_filtered %>% 
  group_by(REEF_NAME, YEAR) %>% 
  summarize(
    site_count = n_distinct(SITE_NO), 
    .groups = "drop"
  )

# Reefs with exactly 3 sites for all 10 years
complete_reefs <- sites_per_year %>% 
  group_by(REEF_NAME) %>%
  summarize(
    year_count = n_distinct(YEAR),
    all_have_3_sites = all(site_count == 3),
    complete_years = year_count == 10
  ) %>%
  filter(all_have_3_sites & complete_years) %>% 
  pull(REEF_NAME)

reef_locations <- data_raw_full %>% 
  distinct(REEF_NAME, LATITUDE, LONGITUDE) %>%
  filter(REEF_NAME %in% complete_reefs) %>%
  group_by(REEF_NAME) %>% 
  slice(1) %>% 
  ungroup()

# Start with northernmost and southernmost reefs
extreme_reefs <- list(
  min_lat_reef = reef_locations$REEF_NAME[which.min(reef_locations$LATITUDE)],
  max_lat_reef = reef_locations$REEF_NAME[which.max(reef_locations$LATITUDE)]
)

selected_reefs <- unique(c(extreme_reefs$min_lat_reef, extreme_reefs$max_lat_reef))
target_reefs <- 15
reefs_needed <- target_reefs - length(selected_reefs)

# Fill in remaining reefs by spreading across longitude
extreme_longitudes <- reef_locations$LONGITUDE[
  reef_locations$REEF_NAME %in% selected_reefs
]
min_long <- min(extreme_longitudes)
max_long <- max(extreme_longitudes)
#skip the first and last 
target_longitudes <- seq(min_long, max_long, length.out = reefs_needed + 2)[
  (1 + 1):(reefs_needed + 1)
]

available_reefs <- reef_locations %>% 
  filter(!REEF_NAME %in% selected_reefs)

additional_reefs <- character(reefs_needed)

for (i in 1:reefs_needed) {
  if(nrow(available_reefs) > 0) {
    distances <- abs(available_reefs$LONGITUDE - target_longitudes[i])
    closest_idx <- which.min(distances)
    closest_reef <- available_reefs$REEF_NAME[closest_idx]
    
    additional_reefs[i] <- closest_reef
    available_reefs <- available_reefs[-closest_idx, ]
  }
}

selected_reefs <- c(selected_reefs, additional_reefs[nzchar(additional_reefs)])

selected_reefs <- head(unique(selected_reefs), target_reefs)

# Map reef names to numbers for Stan
reef_lookup <- reef_locations %>%
  filter(REEF_NAME %in% selected_reefs) %>%
  arrange(LONGITUDE) %>%
  mutate(reef_id = 1:n()) %>%
  select(reef_id, REEF_NAME)

# Step 3: Prepare data for Stan
reef_data <- site_coral %>%
  filter(REEF_NAME %in% selected_reefs, YEAR >= 1995, YEAR <= 2004) %>%
  left_join(reef_lookup, by = "REEF_NAME") %>%
  arrange(reef_id, SITE_NO, YEAR) %>%
  group_by(REEF_NAME, SITE_NO) %>%
  mutate(site_id = cur_group_id()) %>%
  ungroup()

n_sites <- max(reef_data$site_id)
unique_years <- sort(unique(reef_data$YEAR))
years_relative <- unique_years - min(unique_years)
n_years <- length(years_relative)

t0 <- 0
ts_raw <- years_relative
ts <- ts_raw + 1e-9  # Tiny number to avoid ODE solver issues

# Create reef mapping for Stan
n_reefs <- max(reef_lookup$reef_id)
reef_map <- reef_data %>%
  distinct(site_id, reef_id) %>%
  arrange(site_id) %>%
  pull(reef_id)

# Transform to matrix format for Stan
obs_matrix <- reef_data %>%
  select(YEAR, site_id, site_coral) %>%
  mutate(Year_Index = match(YEAR, unique_years)) %>%
  pivot_wider(
    id_cols = Year_Index, 
    names_from = site_id, 
    values_from = site_coral, 
    names_prefix = "site_"
  ) %>%
  arrange(Year_Index) %>% 
  select(-Year_Index) %>% 
  as.matrix()

obs_indicator <- matrix(1, nrow = nrow(obs_matrix), ncol = ncol(obs_matrix))
n_obs_total <- sum(obs_indicator)

# Time points for predictions
t_gd <- seq(from = min(ts_raw), to = max(ts_raw), length.out = 100)
n_t_gd <- length(t_gd)

t_gd <- t_gd + 1e-9

stan_data <- list(
  N_years = n_years,
  N_reefs = n_reefs,
  N_sites = n_sites,
  N_obs_total = n_obs_total,
  Y_obs = obs_matrix,
  Y_obs_present = obs_indicator,
  t0 = t0,
  ts = ts,
  reef_map = reef_map,
  rtol = 1e-6,
  atol = 1e-6,
  max_num_steps = 100000,
  no_t_gd = n_t_gd,
  t_gd = t_gd
)

saveRDS(stan_data, data_output_file)

# save reef mapping for plotting
reef_map_for_plotting <- reef_locations %>%
  filter(REEF_NAME %in% selected_reefs) %>%
  arrange(LONGITUDE) %>%
  mutate(reef_id = 1:n()) %>%
  select(reef_id, REEF_NAME, LATITUDE, LONGITUDE)

rds_mapping_output_file <- file.path(output_dir, "reef_id_mapping.rds")
saveRDS(reef_map_for_plotting, rds_mapping_output_file)

cat("Compiling Stan model:", model_file, "\n")
compiled_model <- stan_model(model_file)

# Define initialization function for non-centered parameterization
init_fun_pc <- function() {
  list(
    mu_pop = pmax(1e-4, c(1.6, 0.9, 0.2, 0.35) + rnorm(4, 0, 0.1)),
    tau_pop = pmax(1e-4, rep(0.15, 4) + runif(4, 0, 0.05)),
    z_reef = matrix(rnorm(4 * n_reefs, 0, 0.1), nrow = 4),
    tau_site = pmax(1e-4, rep(0.15, 4) + runif(4, 0, 0.05)),
    z_site = matrix(rnorm(4 * n_sites, 0, 0.1), nrow = 4),
    sigma = pmax(1e-4, rep(0.1, n_sites) + runif(n_sites, -0.05, 0.05)),
    y0_init_C = pmax(0.011, pmin(0.989, rep(0.5, n_sites) + runif(n_sites, -0.1, 0.1)))
  )
}

# Step 5: Run MCMC sampling and save results
fit <- sampling(
  compiled_model,
  data = stan_data,
  seed = 123,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  init = init_fun_pc,
  control = list(max_treedepth = 10)
)

saveRDS(fit, fit_output_file)
