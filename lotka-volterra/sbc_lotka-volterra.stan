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

transformed data {
  real<lower = 0> theta_sim[4];
  theta_sim[1] = normal_rng(1, 0.5);
  theta_sim[3] = normal_rng(1, 0.5);
  theta_sim[2] = normal_rng(0.05, 0.05);
  theta_sim[4] = normal_rng(0.05, 0.05);

	real<lower = 0> sigma_sim[2] ;
	sigma_sim[1] = lognormal_rng(-1, 1);
	sigma_sim[2] = lognormal_rng(-1, 1);
	
	real z_init_sim[2];
	z_init_sim[1] = lognormal_rng(log(10), 1);
	z_init_sim[2] = lognormal_rng(log(10), 1);
	
  int <lower =0> N = 20;
  real ts[N] = to_array_1d([ 1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20.]);
  real z_sim[N, 2] 
    = integrate_ode_rk45(dz_dt, z_init_sim, 0, ts, theta_sim, 
                        rep_array(0.0, 0), rep_array(0, 0), 1e-5, 1e-3, 5e2);
  real y_init_sim[2];
  real<lower = 0> y_sim[N, 2];
  for (k in 1:2) {
    y_init_sim[k] = lognormal_rng(log(z_init_sim[k]), sigma_sim[k]);
    y_sim[ , k] = lognormal_rng(log(z_sim[, k]), sigma_sim[k]);
  }
}

parameters {
	real<lower = 0> theta[4];   // theta = { alpha, beta, gamma, delta }
	real<lower = 0> z_init[2];  // initial population
	real<lower = 0> sigma[2];   // error scale
}

transformed parameters {
  real z[N, 2]
    = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta, 
                        rep_array(0.0, 0), rep_array(0, 0), 1e-5, 1e-3, 5e2);
}

model {
	theta[{1, 3}] ~ normal(1, 0.5);
	theta[{2, 4}] ~ normal(0.05, 0.05);
	sigma ~ lognormal(-1, 1);
	z_init ~ lognormal(log(10), 1);
  for (k in 1:2) {
    y_init_sim[k] ~ lognormal(log(z_init[k]), sigma_sim[k]);
    y_sim[ , k] ~ lognormal(log(z[, k]), sigma_sim[k]);
  }
}

generated quantities {
	int <lower = 0, upper = 1> I_lt_sim[8]  
		= {theta[1] < theta_sim[1], theta[2] < theta_sim[2], theta[3] < theta_sim[3], theta[4] < theta_sim[4], 
		   z_init[1] < z_init_sim[1], z_init[2] < z_init_sim[2], sigma[1] < sigma_sim[1], sigma[2] < sigma_sim[2]};

}