data {
  int<lower=0> N; // total number of observations
  int<lower=0> M; // total number of subjects
  vector[N] yij;  // longitudinal outcome (long format)
  int subj[N]; // Subject ID
  vector<lower=0, upper=1>[N] dij;  // event indicator
  vector<lower=0, upper=1>[N] z_prev; // previous treatment indicator
  vector<lower=0>[N] tij; //measurement time
  vector<lower=0>[N] time2; //time since treatment switch
  vector<lower=0>[N] start; //previous measurement time
  
  int<lower=0>= P; // max number of measurements per individual
  vector<lower=0>[N] period;  // measurement number of this individual
}

transformed data {
  matrix[N, 4] X; //design matrix for fixed effect;
  matrix[N, 4] Z; //design matrix for random effect (location part);
  X[, 1] = rep_vector(1, N);
  X[, 2] = start;
  X[, 3] = z_prev;
  X[, 4] = time2;
  Z = X;
}

parameters {
  real<lower=0, upper=10> shape;
  vector[3] alpha;
  vector[4] beta;
  cholesky_factor_corr[5] L_Omega; //prior correlation
  vector<lower=0>[5] tau;
  matrix[5, M] bz;
  real<lower=0> sigma0;
  real gamma;
  // real<lower=0> sb;
}

transformed parameters {
  vector<lower=0>[N] sigma_ij;
  vector [N] mu_ij;
  matrix[M, 5] bi;
  vector[N] log_haz;
  vector<lower=0>[N] cum_haz;
  vector[N] linear_pred;

  bi = (diag_pre_multiply(tau, L_Omega) * bz)';
  for (k in 1:N) {
      sigma_ij[k] = sigma0 * sqrt(exp(gamma * z_prev[k] + bi[subj[k], 5]));
      // sigma_ij[i] = sigma0;
      mu_ij[k] = X[k] * beta + Z[k] * (bi[subj[k], 1:4])';
      linear_pred[k] = alpha[1] + alpha[2] * mu_ij[k] + alpha[3]*2*log(sigma_ij[k]);
      log_haz[k] = log(shape) + (shape - 1) * log(tij[k]) + linear_pred[k];
      cum_haz[k] = exp(linear_pred[k]) * (pow(tij[k], shape) - pow(start[k], shape));
      
  }
}

model {
//priors
  tau ~cauchy(0,2.5);
  sigma0 ~ normal(0, 10); //page129, stan-ref-2.17.0
  to_vector(bz) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(5);
// log-likelihood
  yij ~ normal(mu_ij, sigma_ij);
  target += sum(dij .* log_haz - cum_haz);
}
