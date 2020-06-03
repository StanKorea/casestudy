functions {
  real[] dz_dt(real t,       // time
               real[] z,     // system state {prey, predator}
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {
    real u = z[1];
    real v = z[2];

    real alpha = theta[1];
    real beta = theta[2];
    real gamma = theta[3];
    real delta = theta[4];

    real du_dt = (alpha - beta * v) * u;
    real dv_dt = (-gamma + delta * u) * v;
    return { du_dt, dv_dt };
  }
}

data {
  int<lower = 0> N;           // num measurements
  real ts[N];                 // measurement times > 0
}



generated quantities {
  real<lower = 0> theta[4];   // theta = { alpha, beta, gamma, delta }
  real<lower = 0> z_init[2];  // initial population
  real<lower = 0> sigma[2];   // error scale
  real z[N, 2]
    = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta,
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-5, 1e-3, 5e2);
  
  // data to generate:
  real y_init[2];             // initial measured population
  real<lower = 0> y[N, 2];    // measured population at measurement times
  
  theta[1] = normal_rng(1, 0.5);
  theta[3] = normal_rng(1, 0.5);
  theta[2] = normal_rng(0.05, 0.05);
  theta[4] = normal_rng(0.05, 0.05);
  sigma[1] = lognormal_rng(-1, 1);
  sigma[2] = lognormal_rng(-1, 1);
  z_init[1] = lognormal_rng(log(10), 1);
  z_init[2] = lognormal_rng(log(10), 1);
  
  for (k in 1:2) {
    y_init[k] = lognormal_rng(log(z_init[k]), sigma[k]);
    y[ , k] = lognormal_rng(log(z[, k]), sigma[k]);
    
  }
}