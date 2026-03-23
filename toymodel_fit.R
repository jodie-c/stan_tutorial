# ============================================================
# Load required packages and set random seed
# ============================================================
library(rstan)
library(posterior)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(123)

# ============================================================
# Generate fake data (2 competing species with Gaussian noise)
# ============================================================

# True parameter values (matching the Stan parameterization:
# p[1] = r1, p[2] = r2, p[3] = K)
r1_true    <- 0.15
r2_true    <- 0.10
K_true     <- 4.0
y0_true    <- c(1.2, 0.8)   # initial abundances for the two species
sigma_true <- 0.2           # observation noise standard deviation

# Experimental design: replicate wells and observation times
no_wells <- 6
t0 <- 0
ts <- seq(0, 48, by = 6)
if (ts[1] <= t0) ts[1] <- t0 + 1e-8
no_ts <- length(ts)

# Dense time grid used for posterior predictive visualization
ts_gen <- seq(min(ts), max(ts), length.out = 200)
no_ts_gen <- length(ts_gen)

# Deterministic (noise-free) trajectories for the two-species competition model.
# There is no simple closed-form solution when r1 != r2, so we simulate the
# coupled ODE system numerically using a fourth-order Runge-Kutta integrator.

compete_rhs <- function(t, y, r1, r2, K) {
  y1 <- y[1]
  y2 <- y[2]
  
  c(
    r1 * y1 * (1 - (y1 + y2) / K),
    r2 * y2 * (1 - (y1 + y2) / K)
  )
}

rk4_sim <- function(times, y0, r1, r2, K) {
  y <- matrix(NA_real_, nrow = length(times), ncol = 2)
  y[1, ] <- y0
  
  for (i in 2:length(times)) {
    dt <- times[i] - times[i - 1]
    t  <- times[i - 1]
    yi <- y[i - 1, ]
    
    k1 <- compete_rhs(t,          yi,             r1, r2, K)
    k2 <- compete_rhs(t + 0.5*dt, yi + 0.5*dt*k1, r1, r2, K)
    k3 <- compete_rhs(t + 0.5*dt, yi + 0.5*dt*k2, r1, r2, K)
    k4 <- compete_rhs(t + dt,     yi + dt*k3,     r1, r2, K)
    
    y[i, ] <- yi + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
    y[i, ] <- pmax(y[i, ], 0)  # enforce nonnegative abundances
  }
  y
}

# Noise-free trajectories evaluated at the observation times (dimension: no_ts x 2)
mu_true_ts <- rk4_sim(ts, y0_true, r1_true, r2_true, K_true)

# Plot noise-free trajectories
plot(ts, mu_true_ts[, 1], type = "l", lwd = 2,
     xlab = "time", ylab = "population",
     ylim = c(0, 4.5), xlim = c(0, 48),
     main = "Simulated 2-species competition data")
lines(ts, mu_true_ts[, 2], lwd = 2)

# Generate replicate observations: each well is (mu + N(0, sigma)) for both species
# Stan expects: array[no_wells] matrix[2, no_ts]
y_tilde_list <- vector("list", no_wells)

for (w in 1:no_wells) {
  y1_obs <- rnorm(no_ts, mean = mu_true_ts[, 1], sd = sigma_true)
  y2_obs <- rnorm(no_ts, mean = mu_true_ts[, 2], sd = sigma_true)
  y_tilde_list[[w]] <- rbind(y1_obs, y2_obs)  # 2 x no_ts
  
  points(ts, y1_obs, pch = 16)
  points(ts, y2_obs, pch = 16)
}

# ============================================================
# Run Bayesian inference with Stan
# ============================================================

stan_file <- "~/toymodel.stan"

stan_data <- list(
  no_ts     = no_ts,
  no_wells  = no_wells,
  t0        = t0,
  ts        = ts,
  y_tilde   = y_tilde_list,
  no_ts_gen = no_ts_gen,
  ts_gen    = ts_gen
)

fit <- stan(
  file = stan_file,
  data = stan_data,
  chains = 4,
  iter = 2000,
  seed = 123
)

# ============================================================
# Extract draws from fit
# ============================================================
post <- rstan::extract(fit)

# ============================================================
# Print summary of fit
# ============================================================
print(fit, pars = c("p", "y0", "sigma"), probs = c(0.05, 0.5, 0.95))

# ============================================================
# Inspect pair plots
# ============================================================
mcmc_pairs(fit, pars = c("p[1]", "p[2]", "p[3]", "y0[1]", "y0[2]", "sigma"), 
           off_diag_args = list(size = 0.5, alpha = 0.2), 
           diag_fun = "dens", diag_args = list(alpha = 0.6))
