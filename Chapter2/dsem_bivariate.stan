
  data {
  int<lower = 1> N;
  int<lower = 1> N_subj;
  array[N] int<lower=1, upper=N_subj> subj;
  vector[N] YA;
  vector[N] YB;
  vector[N] time;
  vector[N] YA_tm1;
  vector[N] YB_tm1;
  vector[N_subj] mtime;
}

parameters {
  real sigma1;
  real sigma2;
  vector<lower=0>[10] tau_u;
  real alpha1;
  real beta1;
  real phi1;
  real alpha2;
  real beta2;
  real phi2;
  real couple1;
  real couple3;
  matrix[10, N_subj] z_u;
  cholesky_factor_corr[10] L_u;
}

transformed parameters {
matrix[N_subj,10] u;
u = transpose(diag_pre_multiply(tau_u, L_u) * z_u);
}

model {
  // priors EEG
  alpha1 ~ normal(0, 50);
  beta1 ~ normal(0,50);
  sigma1 ~ normal(0, 50);
  phi1 ~ normal(0,50);
  couple1 ~ normal(0,50);
  // priors RT
  alpha2 ~ normal(0, 50);
  beta2 ~ normal(0,50);
  sigma2 ~ normal(0, 50);
  phi2 ~ normal(0,50);
  couple3 ~ normal(0,50);
  tau_u ~ cauchy(0,2);
  to_vector(z_u) ~ std_normal();
  L_u ~ lkj_corr_cholesky(1);
  // Specify model
  YA ~ normal(alpha1 + u[subj,1] + (time - mtime[subj]).*(beta1+u[subj,2]) + (YA_tm1 - (alpha1 + u[subj,1])).*(phi1+u[subj,3]) + (YB_tm1 - (alpha2 + u[subj,6])).*(couple1+u[subj,5]), 
              exp(sigma1 + u[subj,4]));
  YB ~ normal(alpha2 + u[subj,6] + (time - mtime[subj]).*(beta2+u[subj,7]) + (YB_tm1 - (alpha2 + u[subj,6])).*(phi2+u[subj,8]) + (YA_tm1 - (alpha1 + u[subj,1])).*(couple3+u[subj,10]), 
              exp(sigma2 + u[subj,9]));
}

generated quantities {
corr_matrix[10] rho_u = L_u * L_u';

}

