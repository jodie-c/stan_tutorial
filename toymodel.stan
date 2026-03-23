functions {
  vector logistic_ode(real t, vector y, vector p) {
    vector[2] dydt;
    real r1 = p[1];
    real r2 = p[2];
    real K  = p[3];
    
    // ---- ODE model ----
    dydt[1] = r1 * y[1] * (1 - (y[1] + y[2]) / K);
    dydt[2] = r2 * y[2] * (1 - (y[1] + y[2]) / K);

    return dydt;
  }
}

data {
  int<lower=1> no_ts;                           // Number of observation time points
  int<lower=1> no_wells;                        // Number of wells
  real t0;                                      // Initial time
  array[no_ts] real<lower=t0> ts;               // Observation times
  array[no_wells] matrix[2, no_ts] y_tilde;     // Observed data

  int<lower=1> no_ts_gen;                       // Number of prediction time points
  array[no_ts_gen] real ts_gen;                 // Prediction times
}

transformed data {
  //transformed data here
}

parameters {
  vector<lower=0>[3] p;                         // p[1]=r1, p[2]=r2, p[3]=K
  real<lower=0> sigma;                          // Observation noise
  vector<lower=0>[2] y0;                        // Initial values for both species
}

transformed parameters {
  vector[2] y0_state;                           // Initial state vector
  array[no_ts] vector[2] y_state;               // Deterministic ODE solution (solver output)
  matrix[2, no_ts] y;                           // Deterministic trajectories (for likelihood)

  // ---- Solve ODE system at observation times ----
  y0_state = y0;
  y_state = ode_rk45(logistic_ode, y0_state, t0, ts, p);

  // ---- Convert solver output to 2 x no_ts matrix ----
  for (i in 1:no_ts) {
    y[1, i] = y_state[i][1];
    y[2, i] = y_state[i][2];
  }
}

model {
  // ---- Priors ----
  p[1] ~ normal(0.2, 0.1);                      // Growth rate r1
  p[2] ~ normal(0.2, 0.1);                      // Growth rate r2
  p[3] ~ normal(5, 1);                          // Shared carrying capacity K
  sigma ~ normal(0, 0.5);                       // Noise standard deviation
  y0 ~ normal(1, 0.5);                          // Initial values

  // ---- Likelihood (elementwise) ----
  for (j in 1:no_wells) {
    to_vector(y_tilde[j]) ~ normal(to_vector(y), sigma); 
  }
}

generated quantities {
  // ---- Predictions at time points ts_gen ----
  array[no_ts_gen] vector[2] y_state_pred;      // Deterministic ODE solution (solver output)
  matrix[2, no_ts_gen] y_mean;                  // Noise-free trajectory 
  matrix[2, no_ts_gen] y_pred;                  // Posterior predictive draws (with noise)

  vector[2] y0_state_gen;
  y0_state_gen = y0;

  // ---- Solve ODE at prediction times ----
  y_state_pred = ode_rk45(logistic_ode, y0_state_gen, t0, ts_gen, p);

  // ---- Store noise-free mean + generate noisy observations ----
  for (i in 1:no_ts_gen) {
    y_mean[1, i] = y_state_pred[i][1];         // Noise-free ODE solution
    y_mean[2, i] = y_state_pred[i][2];                 

    y_pred[1, i] = y_mean[1, i] + normal_rng(0, sigma); // Noisy prediction
    y_pred[2, i] = y_mean[2, i] + normal_rng(0, sigma); 
  }
}