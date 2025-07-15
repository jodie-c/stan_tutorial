library(rstan)
library(dplyr)
library(tidyr)
library(magrittr)
library(stats) 

data_file <- "/Users/chengu/Documents/master thesis/new_reef_data.csv"
model_file <- "/Users/chengu/Documents/master thesis/final_final_final_code/stan_file_for_pooling_method/partial_pooling_non_centered_withoutforgetting_0703.stan"
output_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data_output_file <- file.path(output_dir, "stan_data_partial_noncentered_sparse.rds")
fit_output_file <- file.path(output_dir, "stan_fit_partial_noncentered_sparse.rds")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(123) 

# Step 1: Load and clean data
data_raw_full <- read.table(data_file, header = TRUE, sep = ";", 
                           stringsAsFactors = FALSE, quote = "")

# Select relevant columns and process data types
data_raw <- data_raw_full %>%
  select(REEF_NAME, SAMPLE_DATE, SITE_NO, GROUP_CODE, COVER, LATITUDE, LONGITUDE) %>%
  mutate(
    COVER = as.numeric(gsub("[^0-9.]", "", COVER)),
    SAMPLE_DATE = as.Date(SAMPLE_DATE, format = "%Y-%m-%d"),
    YEAR = as.numeric(format(SAMPLE_DATE, "%Y"))
  ) %>%
  filter(!is.na(YEAR), !is.na(COVER))

# Calculate coral cover per site (proportion)
site_coral <- data_raw %>%
  filter(GROUP_CODE %in% c("Hard Coral")) %>%
  group_by(REEF_NAME, YEAR, SITE_NO) %>%
  summarize(
    site_coral = if(all(!is.na(COVER))) sum(COVER, na.rm = TRUE) / 100 else NA, 
    .groups = "drop"
  ) %>%
  mutate(site_coral = ifelse(is.finite(site_coral), site_coral, NA))

# Step 2: Select 15 reefs with complete data from 1995 to 2004
site_coral_filtered <- site_coral %>% 
  filter(YEAR >= 1995, YEAR <= 2004)

sites_per_year <- site_coral_filtered %>% 
  group_by(REEF_NAME, YEAR) %>% 
  summarize(site_count = n_distinct(SITE_NO), .groups = "drop")

# Find reefs with exactly 3 sites for all 10 years
complete_reefs <- sites_per_year %>% 
  group_by(REEF_NAME) %>%
  summarize(
    year_count = n_distinct(YEAR), 
    all_have_3_sites = all(site_count == 3), 
    complete_years = year_count == 10
  ) %>%
  filter(all_have_3_sites & complete_years) %>% 
  pull(REEF_NAME)

# Get unique reef locations
reef_locations <- data_raw_full %>% 
  distinct(REEF_NAME, LATITUDE, LONGITUDE) %>%
  filter(REEF_NAME %in% complete_reefs) %>%
  group_by(REEF_NAME) %>% 
  slice(1) %>% 
  ungroup()

# Step 3: Select 15 reefs spread across latitude and longitude
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
target_longitudes <- seq(min_long, max_long, length.out = reefs_needed + 2)[
  (1 + 1):(reefs_needed + 1)
]

available_reefs <- reef_locations %>% 
  filter(!REEF_NAME %in% selected_reefs)

additional_reefs <- character(reefs_needed)
available_indices <- seq_len(nrow(available_reefs))

for (i in 1:reefs_needed) {
  distances <- abs(available_reefs$LONGITUDE[available_indices] - target_longitudes[i])
  closest_rel_idx <- which.min(distances)
  closest_abs_idx <- available_indices[closest_rel_idx]
  
  closest_reef <- available_reefs$REEF_NAME[closest_abs_idx]
  additional_reefs[i] <- closest_reef
  available_indices <- available_indices[-closest_rel_idx]
}

selected_reefs <- c(selected_reefs, additional_reefs[nzchar(additional_reefs)])
selected_reefs <- head(selected_reefs, target_reefs)

# Step 4: Map reef names to numbers for Stan
reef_lookup <- reef_locations %>%
  filter(REEF_NAME %in% selected_reefs) %>%
  arrange(LONGITUDE) %>%
  mutate(Reef_ID = 1:n()) %>%
  select(Reef_ID, REEF_NAME)

# Step 5: Prepare initial data with all sites
all_sites_data <- site_coral %>%
  filter(REEF_NAME %in% selected_reefs, YEAR >= 1995, YEAR <= 2004) %>%
  left_join(reef_lookup, by = "REEF_NAME") %>%
  mutate(SITE_NO = as.character(SITE_NO)) %>%
  arrange(Reef_ID, SITE_NO, YEAR) %>%
  mutate(initial_site_id = match(paste(Reef_ID, SITE_NO, sep="_"), 
                                unique(paste(Reef_ID, SITE_NO, sep="_"))))

# Step 6: Create sparse data structure (5x3, 5x2, 5x1 sites per reef)
# Randomly assign reefs to have 3, 2, or 1 site(s)
unique_reefs <- reef_lookup$REEF_NAME
reefs_3_sites <- sample(unique_reefs, 5)
remaining_reefs <- setdiff(unique_reefs, reefs_3_sites)
reefs_2_sites <- sample(remaining_reefs, 5)
reefs_1_site <- setdiff(remaining_reefs, reefs_2_sites)

# Filter the data to create sparse structure
reduced_data <- all_sites_data %>%
  filter(
    # Keep all sites for reefs_3_sites 
    # Keep only sites '1' and '2' for reefs_2_sites  
    !(REEF_NAME %in% reefs_2_sites & SITE_NO == "3") &
    !(REEF_NAME %in% reefs_1_site & (SITE_NO == "2" | SITE_NO == "3"))
  )

# Step 7: Prepare final data with updated site ids
site_mapping <- reduced_data %>%
  distinct(Reef_ID, SITE_NO) %>%
  arrange(Reef_ID, SITE_NO) %>%
  mutate(final_site_id = 1:n())

final_data <- reduced_data %>%
  left_join(site_mapping, by = c("Reef_ID", "SITE_NO"))

# Recalculate dimensions based on final data
n_sites <- max(final_data$final_site_id)
n_reefs <- max(final_data$Reef_ID)
unique_years <- sort(unique(final_data$YEAR))
years_relative <- unique_years - min(unique_years)
n_years <- length(years_relative)

t0 <- 0
ts <- years_relative + 1e-9  # tiny number to avoid ode solver issues

# Create final reef map for Stan
reef_map <- site_mapping %>%
  arrange(final_site_id) %>%
  pull(Reef_ID)

# Step 8: Transform to matrix format for Stan (with missing data handling)
coral_matrix <- final_data %>%
  select(YEAR, final_site_id, site_coral) %>%
  mutate(Year_Index = match(YEAR, unique_years)) %>%
  pivot_wider(
    id_cols = Year_Index, 
    names_from = final_site_id, 
    values_from = site_coral, 
    names_prefix = "site_"
  ) %>%
  arrange(Year_Index) %>% 
  select(-Year_Index) %>% 
  as.matrix()

# Create indicator matrix for missing data
obs_indicator <- ifelse(is.na(coral_matrix), 0, 1)
n_obs_total <- sum(obs_indicator)

# Replace NAs with placeholder for Stan
coral_matrix[is.na(coral_matrix)] <- 0

t_gd <- seq(from = min(years_relative), to = max(years_relative), length.out = 100)
n_t_gd <- length(t_gd)
t_gd <- pmax(t_gd + 1e-9, ts[1])

# Step 9: Assemble Stan data list
stan_data <- list(
  N_years = n_years,
  N_reefs = n_reefs,
  N_sites = n_sites,
  N_obs_total = n_obs_total,
  Y_obs = coral_matrix,
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

# Save reef mapping for plotting
mapping_output_file <- file.path(output_dir, "reef_id_mapping_sparse.rds")
saveRDS(
  reef_lookup %>% 
  left_join(reef_locations, by="REEF_NAME") %>%
  select(Reef_ID, REEF_NAME, LATITUDE, LONGITUDE),
  mapping_output_file
)

# Step 10: Compile Stan model
compiled_model <- stan_model(model_file)

# Step 11: Run MCMC sampling and save results
fit <- sampling(
  compiled_model,
  data = stan_data,
  seed = 123,
  chains = 4,
  iter = 1000,
  warmup = 500
)

saveRDS(fit, fit_output_file)

