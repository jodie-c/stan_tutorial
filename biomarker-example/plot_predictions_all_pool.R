rm(list = ls())

# Load necessary libraries
library(rstan)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(bayesplot)
library(cowplot)
library(shinystan)
library(gridExtra)
library(grid)
library(loo)


# Script for plotting predictions with generated data. 
# This script generates one pdf with all 10 responsive patients and one pdf with only one user specified patient

#----------------- USER INPUTS ---------------------------------------------------------


# Parameters to load the wanted stan-fit.
date <- '2025-04-17' 
experiment <- 'first_cycle' # 'all_but_last_cycle' or 'first_cycle'
number_of_patients <- 10 


# specific patients to plot, patients 1, 4, 13, 14, 15, 17, 24, 26, 28, 29 are available

patient_IDs_to_plot <- c(14, 17, 29)#c(1, 4, 13, 14, 15, 17, 24, 26, 28, 29)

#----------------------------------------------------------------------------------------------

######## create data_in ############
inputdatafile_PSA <- sprintf('./datasets/%s_%d_PSA_df.csv', experiment, number_of_patients)
inputdatafile_day <- sprintf('./datasets/%s_%d_day_df.csv', experiment, number_of_patients)
inputdatafile_Tx <- sprintf('./datasets/%s_%d_Tx_df.csv', experiment, number_of_patients)

# Read data from files:
data_PSA <- read.table(inputdatafile_PSA, header = TRUE, sep = ",", stringsAsFactors = FALSE)
data_day <- read.table(inputdatafile_day, header = TRUE, sep = ",", stringsAsFactors = FALSE)
data_Tx <- read.table(inputdatafile_Tx, header = TRUE, sep = ",", stringsAsFactors = FALSE)

patient_IDs <- c(1, 4, 13, 14, 15, 17, 24, 26, 28, 29, 30, 36, 37, 39, 44, 50, 55, 58, 60, 61)


# Get the number of sampling times & experimental replicates:
no_reps = ncol(data_PSA);
no_ts_max = nrow(data_PSA)-1; #exclude initial row
ts_lengths = colSums(!is.na(data_day %>% slice(-1)))

#remove NA data
data_day[is.na(data_day)] <- 0
data_PSA[is.na(data_PSA)] <- 0
data_Tx[is.na(data_Tx)] <- 0

# Get initial & sampling times from data:
t0_data = data_day %>% slice(1)
ts_data = data_day %>% slice(-1)

t0_data <- as.numeric(unlist(t0_data))

# Get initial & sampling time population-sizes from data:
y1_data_0 = data_PSA %>% slice(1)
y1_data = data_PSA %>% slice(-1)

y1_data_0 <- as.numeric(unlist(y1_data_0))

Tx_data = data_Tx

# for plotting whole series

inputdatafile_PSA_all <- sprintf('./datasets/%d_PSA_df.csv', number_of_patients)
inputdatafile_day_all <- sprintf('./datasets/%d_day_df.csv', number_of_patients)
inputdatafile_Tx_all <- sprintf('./datasets/%d_Tx_df.csv', number_of_patients)

data_PSA_all <- read.table(inputdatafile_PSA_all, header = TRUE, sep = ",", stringsAsFactors = FALSE)
data_day_all <- read.table(inputdatafile_day_all, header = TRUE, sep = ",", stringsAsFactors = FALSE)
data_Tx_all <- read.table(inputdatafile_Tx_all, header = TRUE, sep = ",", stringsAsFactors = FALSE)


# Get the number of sampling times & experimental replicates
no_ts_max_all = nrow(data_PSA_all)-1; #exclude initial row
ts_lengths_all = colSums(!is.na(data_day_all %>% slice(-1)))

#keep only no_reps columns
data_day_all = data_day_all[1:no_reps]
data_PSA_all = data_PSA_all[1:no_reps]
data_Tx_all = data_Tx_all[1:no_reps]

#remove NA data
data_day_all[is.na(data_day_all)] <- 0
data_PSA_all[is.na(data_PSA_all)] <- 0
data_Tx_all[is.na(data_Tx_all)] <- 0


# Get initial & sampling times from data
ts_data_all = data_day_all %>% slice(-1)

Tx_data_all = data_Tx_all

#patient to plot
patient_col <- 1 

patient_ID <- patient_IDs[patient_col]

# Put data in df for plotting
plot_df<- data.frame(
  x = data_day_all[0:ts_lengths_all[patient_col]+1,patient_col],
  y = data_PSA_all[0:ts_lengths_all[patient_col]+1,patient_col]
)


nr_plot_points = nrow(plot_df)-1 # dont count 0

max_days <- max(data_day_all[, patient_col])
t_gd <- 1:max_days

pool_asgmt <- 1:no_reps


# Data to send to Stan
data_in = list(
  no_reps = no_reps,
  no_ts_max = no_ts_max,
  
  t0_data = t0_data,
  ts_lengths = ts_lengths,
  ts_data = transpose(ts_data),
  
  y1_data_0 = y1_data_0,
  y1_data = transpose(y1_data),
  Tx_data = transpose(Tx_data),
  
  #for plotting
  no_ts_max_all = no_ts_max_all, 
  ts_data_all = transpose(ts_data_all),
  Tx_data_all = transpose(Tx_data_all),
  
  no_t_gd = length(t_gd),
  t_gd=t_gd,
  
  patient_col=patient_col,
  nr_plot_points = nr_plot_points
  
)

#read fit

fit_comp_pool <- readRDS(sprintf("./results/fit_%s_%d_comp_pool_%s.rds",experiment, number_of_patients, date))
mod_comp_pool <- stan_model("./models/model_comp_pool_gq.stan")
posterior_draws_comp_pool <- as.matrix(fit_comp_pool, pars = fit_comp_pool@model_pars)
gq_fit_comp_pool <- gqs(mod_comp_pool, draws = posterior_draws_comp_pool, data = data_in)


fit_no_pool <- readRDS(sprintf("./results/fit_%s_%d_no_pool_%s.rds",experiment, number_of_patients, date))
mod_no_pool <- stan_model("./models/model_no_pool_gq.stan")
posterior_draws_no_pool <- as.matrix(fit_no_pool, pars = fit_no_pool@model_pars)
gq_fit_no_pool <- gqs(mod_no_pool, draws = posterior_draws_no_pool, data = data_in)


fit_part_pool <- readRDS(sprintf("./results/fit_%s_%d_part_pool_%s.rds",experiment, number_of_patients, date))
mod_part_pool <- stan_model("./models/model_part_pool_gq.stan")
posterior_draws_part_pool <- as.matrix(fit_part_pool, pars = fit_part_pool@model_pars)
gq_fit_part_pool <- gqs(mod_part_pool, draws = posterior_draws_part_pool, data = data_in)

########## Generate predictive plots ########## 

plots <- list()
i <- 1
  
for(patient_ID_to_plot in patient_IDs_to_plot) {
  
  # plot specific patient 
  patient_col <- match(patient_ID_to_plot, patient_IDs)
  
  patient_ID <- patient_ID_to_plot
  
  ts <- seq(1,ts_data_all[ts_lengths_all[patient_col], patient_col],20) 
  
  # Get Stan-generated predictions complete pooling
  pred_comp_pool <- as.data.frame(gq_fit_comp_pool,
                         pars = grep(sprintf("y_fit_tot\\[.*,%d\\]", patient_col),
                                     gq_fit_comp_pool@sim$fnames_oi,
                                     value = TRUE)) %>%
    gather(factor_key = TRUE) %>%  # Reshape data frame
    group_by(key) %>%
    summarize(lb = quantile(value, probs = 0.025),
              lb_25 = quantile(value, probs = 0.25),
              median = quantile(value, probs = 0.5),
              ub_75 = quantile(value, probs = 0.75),
              ub = quantile(value, probs = 0.975))
  
  pred_comp_pool_subset <- pred_comp_pool[ts,]
  pred_comp_pool_subset$key <- ts
  
  
  # Get Stan-generated predictions no pooling
  pred_no_pool <- as.data.frame(gq_fit_no_pool,
                                pars = grep(sprintf("y_fit_tot\\[.*,%d\\]", patient_col),
                                            gq_fit_no_pool@sim$fnames_oi,
                                            value = TRUE)) %>%
    gather(factor_key = TRUE) %>%  # Reshape data frame
    group_by(key) %>%
    summarize(lb = quantile(value, probs = 0.025),
              lb_25 = quantile(value, probs = 0.25),
              median = quantile(value, probs = 0.5),
              ub_75 = quantile(value, probs = 0.75),
              ub = quantile(value, probs = 0.975))
  
  pred_no_pool_subset <- pred_no_pool[ts,]
  pred_no_pool_subset$key <- ts
  
  # Get Stan-generated predictions partial pooling
  pred_part_pool <- as.data.frame(gq_fit_part_pool,
                         pars = grep(sprintf("y_fit_tot\\[.*,%d\\]", patient_col),
                                     gq_fit_part_pool@sim$fnames_oi,
                                     value = TRUE)) %>%
    gather(factor_key = TRUE) %>%  # Reshape data frame
    group_by(key) %>%
    summarize(lb = quantile(value, probs = 0.025),
              lb_25 = quantile(value, probs = 0.25),
              median = quantile(value, probs = 0.5),
              ub_75 = quantile(value, probs = 0.75),
              ub = quantile(value, probs = 0.975))
  
  pred_part_pool_subset <- pred_part_pool[ts,]
  pred_part_pool_subset$key <- ts
  
  # Plot prediction intervals
  tsplot_cp <- ggplot(pred_comp_pool_subset, aes(x = as.numeric(key), y = median)) +
    
    geom_vline(xintercept = ts_data[ts_lengths[patient_col], patient_col],
               linetype = "dashed", color = "black", linewidth = 0.5) +
  
    # # complete pooling prediction
    # geom_ribbon(aes(ymin = lb, ymax = ub), fill = "chartreuse3", alpha = 0.1) +
    # geom_line(color = "chartreuse3", linewidth = 1) +
    # geom_line(aes(y = lb), color = "chartreuse3", alpha = 0.5) +
    # geom_line(aes(y = ub), color = "chartreuse3", alpha = 0.5) +
  
    # No pooling prediction
    geom_ribbon(data = pred_no_pool_subset,
                aes(x = as.numeric(key), ymin = lb, ymax = ub),
                fill = "deepskyblue4", alpha = 0.1,) +
    geom_line(data = pred_no_pool_subset,
              aes(x = as.numeric(key), y = median),
              color = "deepskyblue4", linewidth = 1) +
    geom_line(data = pred_no_pool_subset,
              aes(x = as.numeric(key), y = lb),
              color = "deepskyblue4", alpha = 0.5) +
    geom_line(data = pred_no_pool_subset,
              aes(x = as.numeric(key), y = ub),
              color = "deepskyblue4", alpha = 0.5) +
    
    # part pooling prediction
    geom_ribbon(data = pred_part_pool_subset,
                aes(x = as.numeric(key), ymin = lb, ymax = ub),
                fill = "brown2", alpha = 0.1) +
    geom_line(data = pred_part_pool_subset,
              aes(x = as.numeric(key), y = median),
              color = "brown2", linewidth = 1) +
    geom_line(data = pred_part_pool_subset,
              aes(x = as.numeric(key), y = lb),
              color = "brown2", alpha = 0.5) +
    geom_line(data = pred_part_pool_subset,
              aes(x = as.numeric(key), y = ub),
              color = "brown2", alpha = 0.5) + 
    coord_cartesian(ylim = c(0, 30), xlim=c(0,ts[length(ts)]),expand=FALSE)
  
  
  # Put data in df for plotting
  plot_df <- data.frame(
    x = data_day_all[0:ts_lengths_all[patient_col] + 1, patient_col],
    y = data_PSA_all[0:ts_lengths_all[patient_col] + 1, patient_col]
  )
  
  
  # Plot data as points
  tsplot_cp <- tsplot_cp +
    geom_point(data = plot_df, aes(x = x, y = y), show.legend = FALSE, size = 1.2, color = "grey50") 
  
  # Add labels
  tsplot_cp <- tsplot_cp +
    # labs(y = "PSA (ug/L)",x = "Time (days)") +
    annotate("text",x = -Inf, y = Inf, hjust = -0.75, vjust = 3, size = 6,label=paste0("Patient ",patient_ID),
              family="Times", fontface = "bold") +
    theme_bw(base_size = 16) +
    theme(text=element_text(family="Times"),
          panel.grid.major = element_line(colour = "grey90"),
          axis.title = element_blank(),
          panel.grid.minor = element_blank())
  
  plots[[i]] <- tsplot_cp
  i <- i + 1

}

# plots[[1]] <- plots[[1]] +
#   annotate("text",x = -Inf, y = Inf, hjust = -0.75, vjust = 1.5, size = 6,label="(b)",
#            family="Times", fontface = "bold")

# pdf(sprintf("./results/output_%s_%d_all_pool_%s.pdf", experiment, number_of_patients, date),  width=7, height=15)
pdf(sprintf("./results/output_%s_%d_no_part_pool_%s.pdf", experiment, number_of_patients, date),  width=7, height=10)

grid.arrange(grobs = plots, ncol = 1,
             bottom = textGrob("Time (days)",
                               gp = gpar(fontsize = 16, fontfamily="Times")),
             left = textGrob("PSA (ug/L)",
                             rot = 90,
                             gp = gpar(fontsize = 16, fontfamily="Times")))

dev.off()


### ELPD-loo calculation 
y_fit_draws_no_pool_llh <- as.array(gq_fit_no_pool, pars = "log_lik") 
y_fit_draws_comp_pool_llh <- as.array(gq_fit_comp_pool, pars = "log_lik") 
y_fit_draws_part_pool_llh <- as.array(gq_fit_part_pool, pars = "log_lik") 

for(patient_col in c(1:10)){
  print(patient_col)
  
  y_fit_draws_subset_no_pool   <- y_fit_draws_no_pool_llh[,,(sum(ts_lengths[1:patient_col-1])+1):sum(ts_lengths[1:patient_col]) ]
  y_fit_draws_subset_part_pool <- y_fit_draws_part_pool_llh[,,(sum(ts_lengths[1:patient_col-1])+1):sum(ts_lengths[1:patient_col]) ]
  y_fit_draws_subset_comp_pool <- y_fit_draws_comp_pool_llh[,,(sum(ts_lengths[1:patient_col-1])+1):sum(ts_lengths[1:patient_col]) ]
  
  y_fit_draws_subset_no_pool   <- y_fit_draws_no_pool_llh[,,(sum(ts_lengths[1:patient_col-1])+1):sum(ts_lengths[1:patient_col]) ]
  y_fit_draws_subset_part_pool <- y_fit_draws_part_pool_llh[,,(sum(ts_lengths[1:patient_col-1])+1):sum(ts_lengths[1:patient_col]) ]
  y_fit_draws_subset_comp_pool <- y_fit_draws_comp_pool_llh[,,(sum(ts_lengths[1:patient_col-1])+1):sum(ts_lengths[1:patient_col]) ]
  
  loo_estimate_test_no_pool   <- loo(y_fit_draws_subset_no_pool)
  loo_estimate_test_comp_pool <- loo(y_fit_draws_subset_comp_pool)
  loo_estimate_test_part_pool <- loo(y_fit_draws_subset_part_pool)
  
  print(loo_estimate_test_no_pool)
  print(loo_estimate_test_no_pool$diagnostics$pareto_k[loo_estimate_test_no_pool$diagnostics$pareto_k > 0.7])
  print(data_in$y1_data[patient_col,loo_estimate_test_no_pool$diagnostics$pareto_k > 0.7])

  print(loo_estimate_test_comp_pool)
  print(loo_estimate_test_comp_pool$diagnostics$pareto_k[loo_estimate_test_comp_pool$diagnostics$pareto_k > 0.7])
  print(data_in$y1_data[patient_col,loo_estimate_test_comp_pool$diagnostics$pareto_k > 0.7])

  print(loo_estimate_test_part_pool)
  print(loo_estimate_test_part_pool$diagnostics$pareto_k[loo_estimate_test_part_pool$diagnostics$pareto_k > 0.7])
  print(data_in$y1_data[patient_col,loo_estimate_test_part_pool$diagnostics$pareto_k > 0.7])
  
  print(loo_compare(loo_estimate_test_no_pool,loo_estimate_test_comp_pool,loo_estimate_test_part_pool))
  
}


loo_estimate_no_pool_all <- loo(y_fit_draws_no_pool_llh)
loo_estimate_part_pool_all <- loo(y_fit_draws_part_pool_llh)
loo_estimate_comp_pool_all <- loo(y_fit_draws_comp_pool_llh)

print(loo_estimate_no_pool_all)
print(loo_estimate_part_pool_all)
print(loo_estimate_comp_pool_all)

print(loo_compare(loo_estimate_no_pool_all,loo_estimate_part_pool_all,loo_estimate_comp_pool_all))

# Create diagnostic plots
pareto_k_no_pool <- pareto_k_values(loo_estimate_no_pool_all)
pareto_k_part_pool <- pareto_k_values(loo_estimate_part_pool_all)
pareto_k_comp_pool <- pareto_k_values(loo_estimate_comp_pool_all)

plot_diagnostics <- function(pareto_k, ylabel, save_name){
  df <- data.frame(observation = seq_along(pareto_k),
                   k = pareto_k,
                   severity = cut(
                     pareto_k,
                     breaks = c(-Inf, 0.7, 1, Inf),
                     labels = c("Good", "OK", "Bad")))
  
  p <- ggplot(df, aes(x = observation, y = k, alpha = severity)) +
    geom_point(size = 5, shape = "+", color="royalblue4") +
    geom_hline(yintercept = c(0.0, 0.7, 1),
               linetype = c("dotted","dashed","solid"),
               color = c("grey", "sandybrown", "orangered4")) +
    scale_alpha_manual(values = c("Good" = 0.4,"OK" = 0.6,"Bad" = 1.0))+
    labs(x = "Observation", y = ylabel) +
    theme_classic(base_size=16) + 
    ylim(-0.4,1.4) +
    theme(legend.position = "none")
  
  ggsave(save_name, p, width = 7, height = 5)
}

plot_diagnostics(pareto_k_no_pool, ylabel="Pareto k value - No pooling", save_name=paste0("./results/",experiment,"_diagnostics_no_pool.pdf"))
plot_diagnostics(pareto_k_comp_pool, ylabel="Pareto k value - Complete pooling", save_name=paste0("./results/",experiment,"_diagnostics_comp_pool.pdf"))
plot_diagnostics(pareto_k_part_pool, ylabel="Pareto k value - Partial pooling", save_name=paste0("./results/",experiment,"_diagnostics_part_pool.pdf"))


