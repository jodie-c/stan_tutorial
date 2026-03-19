rm(list = ls())

#Load libraries and source files for R functions
library(rstan)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(bayesplot)
library(cowplot)
library(shinystan)

######
# Handle input data
######

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


date <- as.character(Sys.Date())

# ---- Parameters to infer the parameters for the desired experiment ---------

experiment <- 'all_but_last_cycle' # 'all_but_last_cycle' or 'first_cycle'

number_of_patients <- 10 

# --------------------------------------------------------------------------

inputdatafile_PSA <- sprintf('./datasets/%s_%d_PSA_df.csv', experiment, number_of_patients)
inputdatafile_day <- sprintf('./datasets/%s_%d_day_df.csv', experiment, number_of_patients)
inputdatafile_Tx <- sprintf('./datasets/%s_%d_Tx_df.csv', experiment, number_of_patients)

# Read data from files
data_PSA <- read.table(inputdatafile_PSA, header = TRUE, sep = ",", stringsAsFactors = FALSE)
data_day <- read.table(inputdatafile_day, header = TRUE, sep = ",", stringsAsFactors = FALSE)
data_Tx <- read.table(inputdatafile_Tx, header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Get the number of sampling times and experimental replicates
no_reps = ncol(data_PSA);
no_ts_max = nrow(data_PSA)-1; #exclude initial row
ts_lengths = colSums(!is.na(data_day %>% slice(-1)))

#remove NA data
data_day[is.na(data_day)] <- -1
data_PSA[is.na(data_PSA)] <- -1
data_Tx[is.na(data_Tx)] <- -1

# Get initial & sampling times from data
t0_data = data_day %>% slice(1)
ts_data = data_day %>% slice(-1)

t0_data <- as.numeric(unlist(t0_data))

# Get initial & sampling time population-sizes from data
y1_data_0 = data_PSA %>% slice(1)
y1_data = data_PSA %>% slice(-1)

y1_data_0 <- as.numeric(unlist(y1_data_0))

Tx_data = data_Tx

# Transform data for later plotting:
inputdata_y1 <- data.frame(populationsize = unlist(y1_data[,1]),replicate = as.vector(matrix(rep(1,length(ts_data[1])),nrow=length(ts_data[1]),byrow=TRUE)),time = rep(ts_data[1],1))

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

# Plotting stuff, not used when sampling

# patient to plot, (patient column and patient ID doesn't match)
patient_col <- 1

# Put data in df for plotting
plot_df<- data.frame(
  x = data_day_all[0:ts_lengths_all[patient_col]+1,patient_col],
  y = data_PSA_all[0:ts_lengths_all[patient_col]+1,patient_col]
)

nr_plot_points = nrow(plot_df)-1 # dont count row 0
max_days <- max(data_day_all[, patient_col])
t_gd <- 1:max_days

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
  
  # for plotting, not used when sampling
  no_ts_max_all = no_ts_max_all, 
  ts_data_all = transpose(ts_data_all),
  Tx_data_all = transpose(Tx_data_all),
  
  no_t_gd = length(t_gd),
  t_gd=t_gd,
  
  patient_col=patient_col,
  nr_plot_points = nr_plot_points
  
)

######
# Compile and run stan file
######

mod1 <- stan_model("./models/model_no_pool.stan")

# Sample
fit1 <- sampling(object = mod1,
                 data = data_in,
                 seed = 1234,
                 chains = 4,
                 iter = 1000,
                 warmup=500
)

saveRDS(fit1, sprintf("./results/fit_%s_%d_no_pool_%s.rds", experiment, number_of_patients, date))

# Print sampling result
print(fit1)
