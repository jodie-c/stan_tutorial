// coral reef ode model with positive priors and fixed initial condition

functions {
  vector coral_rhs(real t, vector y, vector theta) {
    real C = fmax(0.0, fmin(1.0, y[1]));
    real B = fmax(0.0, fmin(1.0, y[2]));
    real avail = fmax(0.0, 1.0 - (C + B));
    vector[2] dydt;
    dydt[1] = theta[1] * C * avail - theta[2] * C + theta[3] * B;
    dydt[2] = theta[2] * C - theta[3] * B - theta[4] * B;
    return dydt;
  }
}

data {
  int<lower=1> N_years;
  int<lower=1> N_sites;
  matrix[N_years, N_sites] Y_obs;
  int<lower=0, upper=1> Y_obs_present[N_years, N_sites];
  real<lower=0> t0;
  array[N_years] real<lower=t0> ts;
  real<lower=0> rtol;
  real<lower=0> atol;
  int<lower=1> max_num_steps;
  int<lower=1> no_t_gd;
  array[no_t_gd] real t_gd;
  real<lower=0, upper=1> N0;  // first observation, passed from R
}

transformed data {
  // fixed initial condition matching MATLAB: C0 = N0^2, B0 = N0*(1-N0)
  vector[2] y0;
  y0[1] = fmax(1e-6, N0 * N0);
  y0[2] = fmax(1e-6, N0 * (1.0 - N0));
}

parameters {
  vector<lower=0>[4] theta;  // strictly positive: alpha, beta, gamma, mu
  vector<lower=0>[N_sites] sigma;
}

model {
  // weakly informative priors, truncated at 0
  theta[1] ~ normal(2.0, 2.5) T[0, ];
  theta[2] ~ normal(1.5, 2.5) T[0, ];
  theta[3] ~ normal(0.3, 1.5) T[0, ];
  theta[4] ~ normal(0.3, 1.5) T[0, ];
  sigma ~ normal(0, 0.25) T[0, ];

  for (s in 1:N_sites) {
    array[N_years] vector[2] yhat = ode_bdf_tol(coral_rhs, y0, t0, ts, rtol, atol, max_num_steps, theta);
    for (t in 1:N_years)
      if (Y_obs_present[t, s] == 1)
        Y_obs[t, s] ~ normal(fmax(0.0, fmin(1.0, yhat[t, 1] + yhat[t, 2])), sigma[s]);
  }
}

generated quantities {
  array[no_t_gd] real y_fit_C;
  array[no_t_gd] real y_fit_B;
  array[no_t_gd] real y_fit_N;

  array[no_t_gd] vector[2] ygd = ode_bdf_tol(coral_rhs, y0, t0, t_gd, rtol, atol, max_num_steps, theta);
  for (g in 1:no_t_gd) {
    y_fit_C[g] = fmax(0.0, fmin(1.0, ygd[g, 1]));
    y_fit_B[g] = fmax(0.0, fmin(1.0, ygd[g, 2]));
    y_fit_N[g] = fmax(0.0, fmin(1.0, ygd[g, 1] + ygd[g, 2]));
  }
}
