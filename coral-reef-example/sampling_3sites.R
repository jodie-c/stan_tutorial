rm(list = ls())

# Load necessary libraries
library(rstan)
library(dplyr)
library(tidyr)
library(readxl)

stan_path  <- "./models/coral_reef_model.stan"
data_path <- "./datasets/new_reef_data.xlsx"
out_dir    <- "./results/"
# dir.create(out_dir, recursive = TRUE)

san <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
clamp01 <- function(x) pmin(1, pmax(0, x))
toPct <- function(x) pmin(100, pmax(0, 100 * x))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
sm <- stan_model(file = stan_path)

all_data <- read_excel(data_path) %>%
  filter(GROUP_CODE == "Hard Coral") %>%
  mutate(YEAR = as.integer(format(as.Date(SAMPLE_DATE), "%Y")), site_coral = COVER / 100) %>%
  group_by(REEF_NAME, YEAR, SITE_NO) %>%
  summarise(site_coral = mean(site_coral), .groups = "drop")

reef_configs <- list(
  list(name = "North Direction Reef", pattern = "NORTH DIRECTION", data_name = "NORTH DIRECTION REEF",  y0 = 1993, y1 = 2012, sites = 1:3),
  list(name = "Chicken Reef",         pattern = "CHICKEN",         data_name = "CHICKEN REEF",          y0 = 2009, y1 = 2023, sites = 1:3),
  list(name = "Turner Reef",          pattern = "TURNER",          data_name = "TURNER REEF",           y0 = 2002, y1 = 2017, sites = 1:3),
  list(name = "Chinaman Reef",        pattern = "CHINAMAN",        data_name = "CHINAMAN REEF(22102)",  y0 = 1992, y1 = 2004, sites = 1:3),
  list(name = "Snake Reef",           pattern = "SNAKE",           data_name = "SNAKE (22088)",         y0 = 1994, y1 = 2008, sites = 1:3),
  list(name = "Horseshoe Reef",       pattern = "HORSESHOE",       data_name = "HORSESHOE",             y0 = 1999, y1 = 2012, sites = 1:3),
  list(name = "Broomfield Reef",      pattern = "BROOMFIELD",      data_name = "BROOMFIELD REEF",       y0 = 2008, y1 = 2022, sites = 1:3),
  list(name = "One Tree Reef",        pattern = "ONE TREE",        data_name = "ONE TREE REEF",         y0 = 1992, y1 = 2006, sites = 1:3),
  list(name = "Lady Musgrave Reef",   pattern = "LADY MUSGRAVE",   data_name = "LADY MUSGRAVE REEF",    y0 = 1992, y1 = 2007, sites = 1:3)
)

for (cfg in reef_configs) {
  dat <- all_data %>%
    filter(grepl(cfg$pattern, REEF_NAME), YEAR >= cfg$y0, YEAR <= cfg$y1, SITE_NO %in% cfg$sites) %>%
    mutate(SITE_MAP = match(SITE_NO, sort(unique(SITE_NO))))
  
  K <- length(unique(dat$SITE_MAP))
  yrs <- sort(unique(dat$YEAR))
  ts <- as.numeric(yrs - min(yrs) + 1e-9)
  t_gd <- as.numeric(seq(0, max(ts), length.out = 400) + 1e-9)
  
  obs_wide <- dat %>%
    mutate(Year_Index = match(YEAR, yrs)) %>%
    pivot_wider(id_cols = Year_Index, names_from = SITE_MAP, values_from = site_coral, names_prefix = "s") %>%
    arrange(Year_Index) %>% select(-Year_Index)
  
  Y_obs <- as.matrix(obs_wide)
  Y_obs[is.na(Y_obs)] <- 0
  Y_mask <- matrix(as.integer(!is.na(as.matrix(obs_wide))), nrow = nrow(obs_wide))
  
  # get first observation N0 (mean across sites at first year) - matching MATLAB
  first_year_data <- dat %>% filter(YEAR == min(yrs))
  N0 <- mean(first_year_data$site_coral, na.rm = TRUE)
  N0 <- max(0.01, min(0.99, N0))
  
  stan_data <- list(N_years = length(yrs), N_sites = K, Y_obs = Y_obs, Y_obs_present = Y_mask,
                    t0 = 0, ts = ts, rtol = 1e-6, atol = 1e-6, max_num_steps = 100000,
                    no_t_gd = length(t_gd), t_gd = t_gd, N0 = N0)
  
  fit <- sampling(sm, data = stan_data, seed = 123, chains = 4, iter = 2000, warmup = 1000,
                  control = list(adapt_delta = 0.9, max_treedepth = 12))
  
  # Save sampling as RDS files
  saveRDS(stan_data, file.path(out_dir, paste0("data_", san(cfg$name), ".rds")))
  saveRDS(fit, file.path(out_dir, paste0("fit_", san(cfg$name), ".rds")))
  
  # Save summary of sampling
  post <- rstan::extract(fit, pars = c("y_fit_N", "y_fit_C", "y_fit_B"))
  mN <- colMeans(post$y_fit_N); sN <- apply(post$y_fit_N, 2, sd)
  mC <- colMeans(post$y_fit_C); mB <- colMeans(post$y_fit_B)
  
  df_pred <- data.frame(
    t_rel = t_gd, year = t_gd + cfg$y0,
    N_mean = toPct(mN), N_sd = 100 * sN,
    N_lo1 = toPct(clamp01(mN - sN)),   N_hi1 = toPct(clamp01(mN + sN)),
    N_lo2 = toPct(clamp01(mN - 2*sN)), N_hi2 = toPct(clamp01(mN + 2*sN)),
    N_lo3 = toPct(clamp01(mN - 3*sN)), N_hi3 = toPct(clamp01(mN + 3*sN)),
    C_mean = toPct(mC), B_mean = toPct(mB)
  )
  
  obs_data <- all_data %>%
    filter(REEF_NAME == cfg$data_name, YEAR >= cfg$y0, YEAR <= cfg$y1) %>%
    transmute(year = YEAR, site = SITE_NO, cover = toPct(site_coral))
  
  theta <- rstan::extract(fit, "theta")$theta
  param_summary <- data.frame(
    parameter = c("alpha", "beta", "gamma", "mu"),
    mean = colMeans(theta), sd = apply(theta, 2, sd),
    q025 = apply(theta, 2, quantile, 0.025), q50 = apply(theta, 2, quantile, 0.5),
    q975 = apply(theta, 2, quantile, 0.975)
  )
  
  reef_safe <- san(cfg$name)
  write.csv(df_pred, file.path(out_dir, paste0(reef_safe, "_predictions.csv")), row.names = FALSE)
  write.csv(obs_data, file.path(out_dir, paste0(reef_safe, "_observations.csv")), row.names = FALSE)
  write.csv(param_summary, file.path(out_dir, paste0(reef_safe, "_parameters.csv")), row.names = FALSE)
  
}

# Save summary of parameter estimates for all reefs
summary_list <- lapply(reef_configs, function(cfg) {
  fit <- readRDS(file.path(out_dir, paste0("fit_", san(cfg$name), ".rds")))
  theta <- as.data.frame(rstan::extract(fit, "theta")$theta)
  colnames(theta) <- c("alpha", "beta", "gamma", "mu")
  data.frame(Reef = cfg$name, t(colMeans(theta)))
})
write.csv(bind_rows(summary_list), file.path(out_dir, "parameter_summary_all_reefs.csv"), row.names = FALSE)

