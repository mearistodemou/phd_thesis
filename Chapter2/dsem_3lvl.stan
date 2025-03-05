data {
  int<lower = 1> N; //number of observations (rows)
  int<lower = 1> N_subj; //number of subjects
  int<lower = 1> N_days; //number of days (maximum)
  array[N] int<lower=1, upper=N_subj> subj; //column of subj numbers
  array[N] int<lower=1, upper=N_days> day; //column of days
  vector[N] Y; //outcome variable
  vector[N] time; //time variable
  vector[N] Y_tm1d; //lagged outcome variable
  vector[N_subj] mtime;
}

parameters {
  real sigma; //residual variance
  vector<lower=0>[4] tau_u; //sd of subject deviations
  real<lower=0> tau_d; //sd of daily deviations
  real alpha; //mean
  real beta; //trend wp
  //real beta2; // trend bp
  real phi; //autoregression (ar1)
  vector[N_days] z_d;
  matrix[4, N_subj] z_u;
  //cholesky correlation matrix, again to ease sampling
  cholesky_factor_corr[4] L_u;
}

transformed parameters {
matrix[N_subj,4] u;
vector[N_days] d;
u = transpose(diag_pre_multiply(tau_u, L_u) * z_u);
d = z_d * tau_d;
}

model {
  alpha ~ normal(0,50);
  phi ~ normal(0,50);
  beta ~ normal(0,50);
  //beta2 ~ normal(0,50);
  sigma ~ normal(0,50);
  tau_u ~ cauchy(0,2);
  tau_d ~ cauchy(0,2);
  to_vector(z_u) ~ std_normal();
  to_vector(z_d) ~ std_normal();
  L_u ~ lkj_corr_cholesky(1);
  Y ~ normal(alpha + u[subj,1] + d[day] + (time - mtime[subj]).*(beta+u[subj,2]) + (Y_tm1d - (alpha + u[subj,1] + d[day])).*(phi+u[subj,3]), 
              exp(sigma + u[subj,4]));
}

generated quantities {
corr_matrix[4] rho_u = L_u * L_u';
vector[N] log_lik;
for (i in 1:N){
log_lik[i] = normal_lpdf(Y[i] | alpha + u[subj[i],1] + d[day[i]] + (time[i] - mtime[subj[i]]).*(beta+u[subj[i],2]) + (Y_tm1d[i] - (alpha + u[subj[i],1] + d[day[i]])).*(phi+u[subj[i],3]), 
              exp(sigma + u[subj[i],4]));
} 
real sum_ll = sum(log_lik);
}


