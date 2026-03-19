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
  int<lower=0> ts_lengths[no_reps];  

  array[no_reps] vector [no_ts_max] ts_data;
  
  vector<lower=0>[no_reps] y1_data_0; 
  array[no_reps] vector[no_ts_max] y1_data;
  
  array[no_reps] vector[no_ts_max +1] Tx_data; //Tx0 included
  
  // for generated quantities, not used when sampling
  int<lower=1> no_ts_max_all;
  array[no_reps] vector [no_ts_max_all] ts_data_all;
  array[no_reps] vector[no_ts_max_all +1] Tx_data_all;
  int<lower=1> no_t_gd;
  array[no_t_gd] real t_gd;
  int patient_col;
  int nr_plot_points;

}

parameters {
  real<lower=0, upper=0.1*2> p1;
  real<lower=0, upper=1*2> p2;
  real<lower=0, upper=0.1*2> p3;
  real<lower=0, upper=0.001*3> p4;
  real<lower=0, upper=1*2> p5;
  real<lower=0> sigma_1;
  real<lower=0> sigma_2;
}

model {
  p1 ~ normal(0.0278, 0.1); //Ps
  p2 ~ normal(0.69, 1); //lambda
  p3 ~ normal(0.036, 0.1); //alpha
  p4 ~ normal(0.000187,0.001); //rho
  p5 ~ normal(0.0856, 1); //phi

  sigma_1 ~normal(0, 1);
  sigma_2 ~normal(0, 1);
//   
  vector[5] p;
  p[1] = p1;
  p[2] = p2;
  p[3] = p3;
  p[4] = p4;
  p[5] = p5;


  vector[3] SDP0;
  for(x in 1:no_reps){
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
      y1_data[x,t] ~normal(mu[t][3], sigma_1 + mu[t][3]*sigma_2 );
    }
  }
}
