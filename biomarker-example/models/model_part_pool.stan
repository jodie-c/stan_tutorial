functions {
  real lookup(array[] real ts, vector Tx_data, real t, int no_ts) {

    if (t < ts[1]){
        return Tx_data[1];
      }

    for (i in 1:(no_ts - 1)) {
      if (t >= ts[i] && t < ts[i + 1]) {
        return Tx_data[i+1];
        }

    }
    return Tx_data[no_ts+1];  // If t is out of bounds, return last value
  }

  vector pca_model(real t,
                   vector SDP,
                   vector p,
                   vector Tx_data,
                   array[] real ts,
                   int no_ts) {
    vector[3] SDPdt;

    real Tx = lookup(ts, Tx_data, t, no_ts);  // retrieve Tx at time t

    SDPdt[1] = SDP[1] / (SDP[1] + SDP[2]) * p[1] * p[2] * SDP[1];

    SDPdt[2] = (1 - (SDP[1] * p[1]) / (SDP[1] + SDP[2])) * p[2] * SDP[1] - p[3] * SDP[2] * Tx;

    SDPdt[3] = p[4] * SDP[2] - p[5] * SDP[3];

    return SDPdt;
  }
}

data {
  int<lower=1> no_reps;
  int<lower=1> no_ts_max;
  
  vector<lower=0>[no_reps] t0_data;
  int<lower=0> ts_lengths[no_reps];  // Use an integer array instead of a vector

  array[no_reps] vector [no_ts_max] ts_data;
  
  
  vector<lower=0>[no_reps] y1_data_0;
  array[no_reps] vector[no_ts_max] y1_data;
  array[no_reps] vector[no_ts_max +1] Tx_data; //Tx0 included
  
  int<lower=1> no_ts_max_all;
  array[no_reps] vector [no_ts_max_all] ts_data_all;
  array[no_reps] vector[no_ts_max_all +1] Tx_data_all;
  
  int<lower=1> no_t_gd; //generated data
  array[no_t_gd] real t_gd;
  
  int patient_col;
  int nr_plot_points;

}

parameters {

  real <lower=0, upper = 0.1*2> p1_mean;
  real <lower=0> p1_std;
  
  real <lower=0, upper = 1*2> p2_mean;
  real <lower=0> p2_std;

  real <lower=0, upper = 0.1*2> p3_mean;
  real <lower=0> p3_std;

  real <lower=0, upper = 0.001*3> p4_mean;
  real <lower=0> p4_std;

  real <lower=0, upper = 1*2> p5_mean;
  real <lower=0> p5_std;
  
  
  // Individual-level parameters
  array[no_reps] real<lower=0, upper=0.1*2> p1;
  array[no_reps] real<lower=0, upper=1*2> p2;
  array[no_reps] real<lower=0, upper=0.1*2> p3;
  array[no_reps] real<lower=0, upper=0.001*3> p4;
  array[no_reps] real<lower=0, upper=1*2> p5;
  
  real<lower=0> sigma_1;
  real<lower=0> sigma_2;
}


model {
  
  p1_mean ~ normal(0.0278, 0.1);
  p1_std ~ normal(0, 0.1);
  
  p2_mean ~ normal(0.693, 1);
  p2_std ~ normal(0, 1);

  p3_mean ~ normal(0.036, 0.1);
  p3_std ~ normal(0, 0.1);

  p4_mean ~ normal(0.000187, 0.001);
  p4_std ~ normal(0, 0.001);

  p5_mean ~ normal(0.0856, 1);
  p5_std ~ normal(0,1);
  
  sigma_1 ~ normal(0, 1);
  sigma_2 ~ normal(0, 1);
  
  
  for (n in 1:no_reps){
    p1[n] ~ normal(p1_mean, p1_std);
    p2[n] ~ normal(p2_mean, p2_std);
    p3[n] ~ normal(p3_mean, p3_std);
    p4[n] ~ normal(p4_mean, p4_std);
    p5[n] ~ normal(p5_mean, p5_std);
  }
  
  vector[3] SDP0;
  vector[5] p;
  
  
  for(x in 1:no_reps)
  {
    p[1] = p1[x];
    p[2] = p2[x];
    p[3] = p3[x];
    p[4] = p4[x];
    p[5] = p5[x];
    
    
    SDP0[1] = 10;
    SDP0[2] = 1000;
    SDP0[3] = y1_data_0[x];
    
    int ts_slice;
    ts_slice = ts_lengths[x];
    
    // Temporary array to store converted data
    real ts_array[ts_slice];  
    for (ts in 1:ts_slice) {
      ts_array[ts] = ts_data[x, ts];
    }
    
    
    array[ts_slice] vector[3] mu = ode_rk45(pca_model, SDP0, t0_data[x], ts_array, p, Tx_data[x], ts_array, ts_slice);

    for (t in 1:ts_slice){
      y1_data[x,t] ~normal(mu[t][3], sigma_1 + mu[t][3]*sigma_2);
    }
  }
}

