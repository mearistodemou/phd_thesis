data {
  int<lower = 1> N;
  int<lower = 1> N_subj;
  array[N] int<lower=1, upper=N_subj> subj;
  vector[N] Y;
  vector[N] time;
  vector[N] Y_tm1;
  vector[N_subj] mtime;
}

parameters {
  real sigma;
  vector<lower=0>[4] tau_u;
  real alpha;
  real beta;
  real phi;
  matrix[4, N_subj] z_u;
  cholesky_factor_corr[4] L_u;
}

transformed parameters {
matrix[N_subj,4] u;
u = transpose(diag_pre_multiply(tau_u, L_u) * z_u);
}

model {
  alpha ~ normal(0, 50);
  phi ~ normal(0,50);
  beta ~ normal(0,50);
  sigma ~ normal(0, 50);
  tau_u ~ cauchy(0,2);
  to_vector(z_u) ~ std_normal();
  L_u ~ lkj_corr_cholesky(1);
  Y ~ normal(alpha + u[subj,1] + (time - mtime[subj]).*(beta+u[subj,2]) + (Y_tm1 - (alpha + u[subj,1])).*(phi+u[subj,3]), 
              exp(sigma + u[subj,4]));
}

generated quantities {
corr_matrix[4] rho_u = L_u * L_u';
vector[N] log_lik;
for (i in 1:N){
log_lik[i] = normal_lpdf(Y[i] | alpha + u[subj[i],1] + (time[i] - mtime[subj[i]]).*(beta+u[subj[i],2]) + (Y_tm1[i] - (alpha + u[subj[i],1])).*(phi+u[subj[i],3]), 
              exp(sigma + u[subj[i],4]));
} 
real sum_ll = sum(log_lik);
}

