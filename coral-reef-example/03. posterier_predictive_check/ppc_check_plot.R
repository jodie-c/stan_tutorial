library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(viridis) 
library(gridExtra)
library(grid)

data_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/pooling_methods_rds"
data_file <- "/Users/chengu/Documents/master thesis/new_reef_data.csv"
output_dir <- "/Users/chengu/Documents/master thesis/final_final_final_code/plots"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
site_colors <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A") # Red, Blue, Green

# step1: prepare the observed data points for 15 reefs (they are same for all models)
data_raw <- read.table(
  data_file, 
  header = TRUE, 
  sep = ";", 
  stringsAsFactors = FALSE, 
  quote = ""
) 

reef_data <- data_raw %>%
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
site_coral <- reef_data %>%
  # Filter for Hard Coral data only
  filter(GROUP_CODE %in% c("Hard Coral")) %>% 
  # Group by reef, year, and site
  group_by(REEF_NAME, YEAR, SITE_NO) %>%
  # Calculate total coral cover (divide by 100 to get proportion)
  summarize(site_coral = sum(COVER) / 100, .groups = "drop")

site_coral_filtered <- site_coral %>% 
  filter(YEAR >= 1995, YEAR <= 2004)

# Count number of sites per reef per year
sites_per_reef_year <- site_coral_filtered %>% 
  group_by(REEF_NAME, YEAR) %>% 
  summarize(
    site_count = n_distinct(SITE_NO), 
    .groups = "drop"
  )

reefs_with_complete_site_counts <- sites_per_reef_year %>% 
  group_by(REEF_NAME) %>%
  summarize(
    year_count = n_distinct(YEAR),
    all_have_3_sites = all(site_count == 3),
    complete_years = year_count == 10
  ) %>%
  filter(all_have_3_sites & complete_years) %>% 
  pull(REEF_NAME)

reef_locations <- data_raw %>% 
  distinct(REEF_NAME, LATITUDE, LONGITUDE) %>%
  filter(REEF_NAME %in% reefs_with_complete_site_counts) %>%
  group_by(REEF_NAME) %>% 
  slice(1) %>% 
  ungroup()

# Start with northernmost and southernmost reefs
extreme_reefs <- list(
  min_lat_reef = reef_locations$REEF_NAME[which.min(reef_locations$LATITUDE)],
  max_lat_reef = reef_locations$REEF_NAME[which.max(reef_locations$LATITUDE)]
)
selected_reefs_names <- unique(c(extreme_reefs$min_lat_reef, extreme_reefs$max_lat_reef))
target_n_reefs <- 15
reefs_needed <- target_n_reefs - length(selected_reefs_names)

# Fill in remaining reefs by spreading across longitude
if (reefs_needed > 0) {
  lat_extreme_longitudes <- reef_locations$LONGITUDE[
    reef_locations$REEF_NAME %in% selected_reefs_names
  ]
  min_long <- min(lat_extreme_longitudes)
  max_long <- max(lat_extreme_longitudes)
  #skip the first and last 
  target_longitudes <- seq(min_long, max_long, length.out = reefs_needed + 2)[
    (1 + 1):(reefs_needed + 1)
  ]
  
  available_reefs <- reef_locations %>% 
    filter(!REEF_NAME %in% selected_reefs_names)
  
  additional_reefs <- character(reefs_needed)
  for (i in 1:reefs_needed) {
    if(nrow(available_reefs) == 0) break
    distances <- abs(available_reefs$LONGITUDE - target_longitudes[i])
    closest_idx <- which.min(distances)
    closest_reef_name <- available_reefs$REEF_NAME[closest_idx]
    additional_reefs[i] <- closest_reef_name
    available_reefs <- available_reefs[-closest_idx, ]
  }
  selected_reefs_names <- c(selected_reefs_names, additional_reefs[nzchar(additional_reefs)])
}

selected_reefs_names <- head(unique(selected_reefs_names), target_n_reefs)

reef_lookup <- reef_locations %>%
  filter(REEF_NAME %in% selected_reefs_names) %>%
  arrange(LONGITUDE) %>%
  mutate(ReefID = 1:n()) %>%
  select(ReefID, REEF_NAME)

print(reef_lookup)


#step2: load the models and prepare the data for each model
models <- list(
  list(
    name = "Complete Pooling",
    fit_file = file.path(data_dir, "stan_fit_complete_pooling.rds"),
    data_file = file.path(data_dir, "stan_data_complete_pooling.rds")
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

#get the shared data for all models
first_model <- models[[1]]
stan_data_file <- first_model$data_file
stan_data <- readRDS(stan_data_file)
n_sites <- stan_data$N_sites
n_years <- stan_data$N_years
obs_matrix <- stan_data$Y_obs
obs_present <- stan_data$Y_obs_present
ts_obs <- stan_data$ts
t_grid <- stan_data$t_gd
print(paste("Model has", n_sites, "sites."))

#sites are mapped to reefs in a consistent way
reef_map <- ceiling(1:n_sites / 3)  # Sites 1-3 = Reef 1, Sites 4-6 = Reef 2, etc.
print("Created consistent reef mapping (3 sites per reef).")

#step3: prepare the data into long format for plotting
# original Y_obs is a matrix [N_years × N_sites] = [10× 45]:
#      Site0.1  Site0.2  Site0.3  ...  Site0.2
# 0.1   0.3      0.4      0.5     ...   0.6
# 0.2   0.7      0.8      0.9     ...   0.1
# ...   ...      ...      ...     ...   ...
# 0.1   0.2      0.3      0.4     ...   0.5

# final long format result: (450 rows)
# Year_Index | SiteID | Cover | Time
#     0.1    |  0.1   |  0.3  | 0.1  <- Site 0.1, Year 0.1
#     0.2    |  0.1   |  0.7  | 0.2  <- Site 0.1, Year 0.2
#     0.1    |  0.1   |  0.2  | 0.1  <- Site 0.1, Year 0.1
#     0.1    |  0.2   |  0.4  | 0.1  <- Site 0.2, Year 0.1
#     0.2    |  0.2   |  0.8  | 0.2  <- Site 0.2, Year 0.2
#     0.1    |  0.2   |  0.3  | 0.1  <- Site 0.2, Year 0.1



data_long <- data.frame(
  Year_Index = rep(1:n_years, n_sites),
  SiteID = rep(1:n_sites, each = n_years),
  # rep(1:45, each=10) = 1,1,1,...,1, 2,2,2,...,2, 3,3,3,...,3...
  # 111...1, 222...2, 333...3, ... (450 values)
  
  Cover = as.vector(obs_matrix), # column-wise
  
  Time = rep(ts_obs, n_sites) 
  # rep(c(1995,1996,...,2004), 45)
  # 1995,1996,...,2004, 1995,1996,...,2004, ... (450 values)
)


# create a filter to identify and remove missing observations from our dataset(for the sparse data set)
mask_data <- data.frame(
  Year_Index = rep(1:n_years, n_sites),
  SiteID = rep(1:n_sites, each = n_years),
  IsPresent = as.vector(obs_present)
)

observed_data <- data_long %>%
  inner_join(mask_data, by = c("Year_Index", "SiteID")) %>%
  filter(IsPresent == 1) %>%  
  mutate(
    ReefID = reef_map[SiteID],
    SiteWithinReef = ((SiteID - 1) %% 3) + 1
  ) %>%
  left_join(reef_lookup, by = "ReefID") %>%  # Add reef names
  select(Time, Cover, SiteID, ReefID, REEF_NAME, SiteWithinReef)


# step4: generate the plots
plots <- list()

for (model_info in models) {
  
  # Load fit and extract posterior samples
  fit <- readRDS(model_info$fit_file)
  y_fit_samples <- rstan::extract(fit, pars = "y_fit_tot")$y_fit_tot
  
  # Process posterior samples into summary format
  site_summaries <- list()
  for (s in 1:n_sites) {
    site_preds <- y_fit_samples[, s, ]
    site_summary <- apply(site_preds, 2, function(x) {
      quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    })
    site_summaries[[s]] <- data.frame(
      TimeGrid = t_grid,
      SiteID = s,
      ReefID = reef_map[s],
      REEF_NAME = reef_lookup$REEF_NAME[reef_lookup$ReefID == reef_map[s]],
      SiteWithinReef = ((s - 1) %% 3) + 1,
      LowerCI = site_summary[1, ],
      MedianCover = site_summary[2, ],
      UpperCI = site_summary[3, ]
    )
  }
  predicted_data <- bind_rows(site_summaries)
  

  #---------------------------------------------------------------------
  # (this part get help from chatgpt)
  model_plot <- ggplot() +
    # Dashed lines for 95% CrI
    geom_line(data = predicted_data,
              aes(x = TimeGrid, y = LowerCI,
                  group = factor(SiteID), # Group by unique site
                  color = factor(SiteWithinReef)), # Use line color
              linetype = "dashed", linewidth = 0.4, alpha = 0.7) +
    geom_line(data = predicted_data,
              aes(x = TimeGrid, y = UpperCI,
                  group = factor(SiteID), # Group by unique site
                  color = factor(SiteWithinReef)), # Use line color
              linetype = "dashed", linewidth = 0.4, alpha = 0.7) +
    # Solid line for Median Prediction
    geom_line(data = predicted_data,
              aes(x = TimeGrid, y = MedianCover,
                  group = factor(SiteID), # Group by unique site
                  color = factor(SiteWithinReef)), # Use line color
              linewidth = 0.8) + # Make median slightly thicker
    # Observed Data Points (colored by SiteWithinReef)
    geom_point(data = observed_data,
               aes(x = Time, y = Cover,
                   color = factor(SiteWithinReef)), # Use point color
               shape = 16, size = 1.7) + # Solid round shape, slightly larger
    # Facet by actual REEF_NAME, arranged in 3 columns and 5 rows
    facet_wrap(~ REEF_NAME, ncol = 3) + # Use 3 columns for layout
    # Use manual color scales for better distinction
    scale_color_manual(name = "Site within Reef", values = site_colors) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = model_info$name,
      subtitle = "95% credible intervals (dashed) for each site",
      x = "Time (Years)",
      y = "Total Coral Cover"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(size = 12, face = "bold"), # Larger reef names
      axis.title = element_text(size = 14), # Larger axis titles
      axis.text = element_text(size = 12), # Explicit size for axis text
      legend.position = "bottom",
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      legend.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      plot.margin = unit(c(0, 0, 0, 0), "cm"), # Remove plot margins completely
      panel.spacing = unit(0.1, "lines"), # Minimize space between panels
      panel.spacing.x = unit(0.1, "lines"),
      panel.border = element_rect(color = "gray80"),
      aspect.ratio = 0.85  # Slightly wider panels
    )
  
  # Save the individual plot with larger dimensions as PDF for better LaTeX integration
  plot_filename <- file.path(output_dir, paste0("ppc_trajectories_", gsub("[^a-zA-Z0-9_]", "_", model_info$name), ".pdf"))
  ggsave(plot_filename, plot = model_plot, width = 25, height = 15, device = "pdf", limitsize = FALSE) # Saved as PDF
  print(paste("  Plot saved to:", plot_filename))
  
  # Store the plot for the combined figure
  plots[[model_info$name]] <- model_plot
}

# Create a combined plot with all models
if (length(plots) > 0) {
  # Create a title grob
  title_grob <- textGrob(
    "Posterior Predictive Checks Across Pooling Methods",
    gp = gpar(fontface = "bold", fontsize = 20)
  )
  
  # Arrange all plots in a grid with the title
  combined_plot <- grid.arrange(
    grobs = plots,
    ncol = 2,
    top = title_grob,
    bottom = textGrob(
      "Observed data (points) vs. posterior predictions (lines) with 95% credible intervals (dashed)",
      gp = gpar(fontsize = 12)
    )
  )
  
  # Save the combined plot with larger dimensions as PDF for better LaTeX integration
  combined_filename <- file.path(output_dir, "ppc_all_models_combined.pdf")
  ggsave(combined_filename, plot = combined_plot, width = 25, height = 20, device = "pdf", limitsize = FALSE) # Saved as PDF
  print(paste("Combined plot saved to:", combined_filename))
}

print("--- All PPC plots generated with white background styling. ---")
