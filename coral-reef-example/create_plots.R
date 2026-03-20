rm(list = ls())

# Plots comparing bayesian posteriors with frequentist estimates

library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(readr)
library(purrr)
library(latex2exp)

rds_dir    <- "./results/"
data_dir <- "./results"
freq_csv   <- "./results/results_table.csv"
out_dir <- "./results/"

san <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

reef_configs <- list(
  # list(name = "North Direction Reef", pattern = "NORTH DIRECTION", data_name = "NORTH DIRECTION REEF",  y0 = 1993, y1 = 2012, sites = 1:3),
  # list(name = "Chicken Reef",         pattern = "CHICKEN",         data_name = "CHICKEN REEF",          y0 = 2009, y1 = 2023, sites = 1:3),
  # list(name = "Turner Reef",          pattern = "TURNER",          data_name = "TURNER REEF",           y0 = 2002, y1 = 2017, sites = 1:3),
  # list(name = "Chinaman Reef",        pattern = "CHINAMAN",        data_name = "CHINAMAN REEF(22102)",  y0 = 1992, y1 = 2004, sites = 1:3),
  # list(name = "Snake Reef",           pattern = "SNAKE",           data_name = "SNAKE (22088)",         y0 = 1994, y1 = 2008, sites = 1:3),
  # list(name = "Horseshoe Reef",       pattern = "HORSESHOE",       data_name = "HORSESHOE",             y0 = 1999, y1 = 2012, sites = 1:3),
  # list(name = "Broomfield Reef",      pattern = "BROOMFIELD",      data_name = "BROOMFIELD REEF",       y0 = 2008, y1 = 2022, sites = 1:3)
  list(name = "One Tree Reef",        pattern = "ONE TREE",        data_name = "ONE TREE REEF",         y0 = 1992, y1 = 2006, sites = 1:3),
  list(name = "Lady Musgrave Reef",   pattern = "LADY MUSGRAVE",   data_name = "LADY MUSGRAVE REEF",    y0 = 1992, y1 = 2007, sites = 1:3)
)

####### Create parameter plots
param_names <- c("alpha", "beta", "gamma", "mu")
param_labels <- c(alpha = "alpha~-~Intrinsic~Growth", beta = "beta~-~Bleaching", 
                  gamma = "gamma~-~Recovery", mu = "mu~-~Mortality")
    
posteriors <- bind_rows(lapply(reef_configs, function(rf) {
  fit <- readRDS(file.path(rds_dir, paste0("fit_", san(rf$name), ".rds")))
  theta <- as.data.frame(rstan::extract(fit, "theta")$theta)
  colnames(theta) <- param_names
  theta$reef <- rf$name
  theta
}))

freq_results <- read_csv(freq_csv, show_col_types = FALSE)

posteriors_long <- posteriors %>%
  pivot_longer(cols = all_of(param_names), names_to = "parameter", values_to = "value") %>%
  filter(value >= 0) %>%
  mutate(parameter = factor(parameter, levels = param_names, labels = param_labels))

freq_long <- freq_results %>%
  pivot_longer(cols = all_of(param_names), names_to = "parameter", values_to = "mean") %>%
  mutate(lower = case_when(parameter == "alpha" ~ alpha_lo, parameter == "beta" ~ beta_lo,
                           parameter == "gamma" ~ gamma_lo, parameter == "mu" ~ mu_lo),
         upper = case_when(parameter == "alpha" ~ alpha_hi, parameter == "beta" ~ beta_hi,
                           parameter == "gamma" ~ gamma_hi, parameter == "mu" ~ mu_hi),
         reef = case_when(grepl("SNAKE", Reef) ~ "Snake Reef", grepl("HORSESHOE", Reef) ~ "Horseshoe Reef", TRUE ~ Reef),
        parameter = factor(parameter, levels = param_names, labels = param_labels)) %>% 
        filter(reef %in% unique(posteriors_long$reef))

reef_order <- rev(unique(posteriors_long$reef))
posteriors_long$reef <- factor(posteriors_long$reef, levels = reef_order)
freq_long$reef <- factor(freq_long$reef, levels = reef_order)

# combined parameter plot 
label_fig <- data.frame(
  reef = "One Tree Reef",
  parameter = "alpha~-~Intrinsic~Growth",
  lab = "(c)"
)

p <- ggplot(posteriors_long, aes(x = value, y = reef, fill = stat(quantile))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      quantiles = c(0.025, 0.975), scale = 0.85, from = 0, 
                      bandwidth = 0.12, color = "cadetblue4",linewidth = 0.2) +
  scale_fill_manual(name = "Bayesian Posterior", 
                    values = c("paleturquoise3", "cadetblue4", "paleturquoise3")) +
  geom_point(data = freq_long, aes(x = mean, y = reef, shape = "Maximum Likelihood Estimate (MLE)"), 
             color = "coral3", size = 1.8, inherit.aes = FALSE,position = position_nudge(y = -0.03)) +
  geom_errorbarh(data = freq_long, aes(xmin = lower, xmax = upper, y = reef, linetype = "95% Confidence Interval"),
                 color = "coral3", height = 0.15, linewidth = 0.6, inherit.aes = FALSE, position = position_nudge(y = -0.03)) +
  scale_shape_manual(name = "Frequentist", values = c("Maximum Likelihood Estimate (MLE)" = 18)) +
  scale_linetype_manual(name = "Frequentist", values = c("95% Confidence Interval" = "solid")) +
  facet_wrap(~ parameter, scales = "free_x", ncol = 4, labeller = labeller(parameter = label_parsed)) +
  coord_cartesian(xlim = c(0, NA)) +
  labs(x = "Rate of change (per year)") +
  geom_text(data = label_fig,
            aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0, vjust = 2, size = 4, family="Times", fontface = "bold")+
  theme_ridges(font_size = 11) +
  theme(text=element_text(family="Times"),
        axis.title.x = element_text(size = 9, hjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 9, angle = 90,hjust=0),
        strip.text = element_text(size = 10),
        legend.position = "none") +
  scale_y_discrete(expand = c(0, 0)) +
  guides(fill = guide_legend(title.position = "left", nrow = 1, order = 1),
         shape = guide_legend(title.position = "left", order = 2,
                              override.aes = list(color = "coral3", size = 3)),
         linetype = guide_legend(title.position = "left", order = 2,
                                 override.aes = list(color = "coral3")))

ggsave(file.path(out_dir, "coral_reef_param_est.pdf"), p, width = 8, height = 3)
# ggsave(file.path(out_dir, "coral_reef_param_est_others.pdf"), p, width = 8, height = 10)

####### Create prediction plots
## Predictions
pred_bayes <- map_dfr(reef_configs, function(rf) {
  read_csv(file.path(data_dir, paste0(san(rf$name), "_predictions.csv"))) %>%
  mutate(reef_name = rf$name) %>%
  mutate(group = paste0(san(rf$name), "_bayes"))
  })

pred_freq <- map_dfr(reef_configs, function(rf) {
  read_csv(file.path(data_dir, paste0(san(rf$name), "_predictions_freq.csv"))) %>%
  mutate(reef_name = rf$name) %>%
  mutate(group = paste0(san(rf$name), "_freq"))
})

pred_all <- rbind(pred_bayes,pred_freq)

## Observations
obs_bayes <- map_dfr(reef_configs, function(rf) {
  read_csv(file.path(data_dir, paste0(san(rf$name), "_observations.csv"))) %>%
  mutate(reef_name = rf$name) %>%
    mutate(group = paste0(san(rf$name), "_bayes"))
  })

obs_freq <- map_dfr(reef_configs, function(rf) {
  read_csv(file.path(data_dir, paste0(san(rf$name), "_observations.csv"))) %>%
    mutate(reef_name = rf$name) %>%
    mutate(group = paste0(san(rf$name), "_freq"))
})

obs_all <- rbind(obs_bayes,obs_freq)

obs_all$site <- factor(
  obs_all$site,levels = sort(unique(obs_all$site)),
  labels = paste("Site", sort(unique(obs_all$site)))
)

strip_labels <- unique(pred_all[c("reef_name", "group")])

label_df <- data.frame(
  group = unique(pred_all[c("group")]),
  reef_name = unique(pred_all[c("reef_name")])
)

label_fig <- data.frame(
  group = "One_Tree_Reef_freq",
  lab = "(b)"
)

p <- ggplot() +
  
  ## Credible intervals
  geom_ribbon(data = pred_all, aes(year, ymin = N_lo3, ymax = N_hi3),
              fill = "grey92") +
  geom_ribbon(data = pred_all, aes(year, ymin = N_lo2, ymax = N_hi2),
              fill = "grey85") +
  geom_ribbon(data = pred_all, aes(year, ymin = N_lo1, ymax = N_hi1),
              fill = "grey78") +

  ## Mean trajectories
  geom_line(data = pred_all, aes(year, N_mean, colour = "Total N(t)"),
            linewidth = 1) +
  geom_line(data = pred_all, aes(year, C_mean, colour = "Healthy C(t)"),
            linewidth = 0.9, linetype = "dashed") +
  geom_line(data = pred_all, aes(year, B_mean, colour = "Bleaching B(t)"),
            linewidth = 0.9, linetype = "dashed") +
  
  ## Observations
  geom_point(data = obs_all, aes(year, cover, shape = factor(site)),,
    size = 1.5, colour = "grey30") +

  ## Facet by reef
  facet_wrap(~factor(group, levels=sort(unique(pred_all$group),decreasing=TRUE)), ncol=2, scales = "free_x",
             labeller = labeller(group = setNames(strip_labels$reef_name, strip_labels$group))) +
  
  ## Scales
  scale_y_continuous(limits = c(0, 100)) +
  
  scale_colour_manual(name = NULL,
    values = c("Total N(t)"    = rgb(0.3, 0.56, 0.0),
               "Healthy C(t)"  = rgb(0.78, 0.35, 0.35),
               "Bleaching B(t)" = rgb(0.33, 0.6, 0.73),
      setNames(
        scales::hue_pal()(length(unique(obs_all$site))),
        as.character(sort(unique(obs_all$site))))),
    breaks = c("Total N(t)", "Healthy C(t)", "Bleaching B(t)",
      as.character(sort(unique(obs_all$site)))),
    labels = c("Total N(t)", "Healthy C(t)", "Bleaching B(t)",
      paste("Site", as.character(sort(unique(obs_all$site)))))) +
  
  scale_shape_manual(name = NULL,
                     values = c(16, 15, 17)[seq_along(unique(obs_all$site))]) +
  
  ## Labels & theme
  labs(x = "Year", y = "Coral Cover (%)") +
  
  geom_text(data = label_df,
    aes(x = -Inf, y = Inf, label = reef_name),
    hjust = -0.15, vjust = 6, size = 5, family="Times", fontface = "bold")+
  geom_text(data = label_fig,
            aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.5, vjust = 2.5, size = 6, family="Times", fontface = "bold")+

  
  theme_bw(base_size = 16) +
  theme(text=element_text(family="Times"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",#"none",#
        strip.text = element_blank(),
        strip.background = element_blank())


ggsave(file.path(out_dir, "coral_reef_pred.pdf"), p, width = 10, height = 8)
# ggsave(file.path(out_dir, "coral_reef_pred_others.pdf"), p, width = 10, height = 14)

