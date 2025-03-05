data {
  int<lower=1> N; 		 	  // number of subjects
  int<lower=1> T; 		 	  // number of observations
  array[N] vector[T] Y; 	  	  // time series data
  array[N] vector[T] time; 	// time variable
  int<lower=0,upper=N*T> N_miss; 	  // number of missing values
  array[N_miss,2] int ii_miss;		  // position of missing values

  // mean and sd used for imputation
  real y_mean;
  real y_sd;
}

parameters {
  vector[4] gamma; 	      	 	// population-level effects
  array[N] vector[4] u;  	 	// subject-specific deviations

  // for covariance matrix of subject-specific deviations
  cholesky_factor_corr[4] L_Omega; 	// Cholesky factor
  vector<lower=0>[4] tau_Omega;		// vector of standard deviations

  // parameter to model missing data
  array[N_miss] real Y_miss;		
}

transformed parameters {
  // construct covariance matrix from Cholesky factors and standard deviations
  corr_matrix[4] R_Omega; // correlation matrix
  R_Omega = multiply_lower_tri_self_transpose(L_Omega); // R = L* L'
  // quad_form_diag: diag_matrix(tau) * R * diag_matrix(tau)
  cov_matrix[4] Omega = quad_form_diag(R_Omega, tau_Omega);
  
  // subject-specific parameters
  vector[N] mu;  			// mean
  vector[N] phi; 			// autoregression (below mean)
  vector[N] psi; 			// residual variance
  vector[N] trend;			// linear trend

  // object for deviation from mean
  array[N] vector[T] delta;
  array[N] vector[T] delta_time;	// object for deviation from mean time	

  // impute missing values into data matrix
  array[N] vector[T] Y_imp = Y;
  if (N_miss > 0) {
    for (it in 1:N_miss) {
      Y_imp[ii_miss[it][1], ii_miss[it][2]] = Y_miss[it];
    }
  }

  for (i in 1:N) {
    // to obtain subject-specific effects, 
    // sum population-level effects and subjects-specific deviations
    mu[i] = gamma[1] + u[i][1];		// eq.8
    phi[i] = gamma[2] + u[i][2];	// eq.9
        // note gamma[4] and u[i][4] are estimated on a log-scale
    // to assume normal distribution which simplifies estimation
    psi[i] = exp(gamma[4] + u[i][4]);	// eq.11
    trend[i] = gamma[3] + u[i][3];	// 
    // calculate deviation by subtracting mean
    delta[i] = Y_imp[i] - mu[i];	// eq.5
    delta_time[i] = time[i] - trend[i]; // 
  }
}

model {
  // prior distributions
  // .. for the population-level parameters
  target += normal_lpdf(gamma | 0, 3); 
  // .. for the Cholesky factor
  target += lkj_corr_cholesky_lpdf(L_Omega | 1.0);
  // .. for the vector of standard deviations
  target += cauchy_lpdf(tau_Omega | 0, 2.5);
  // .. for imputed data points
  target += normal_lpdf(Y_miss | y_mean, y_sd);

  // likelihood
  for (i in 1:N) {
    // subject-specific deviations
    target += multi_normal_lpdf(u[i] | rep_vector(0,4), Omega);
    
    for (t in 2:T) {
      // data given model parameters
      target += normal_lpdf(delta[i][t] | phi[i] * delta[i][t-1] + trend[i] * delta_time[i][t],
	sqrt(psi[i])); // eq.6
    }
  }
}

generated quantities {
  // obtain posterior predictives (to check model fit)
  array[N] vector[T] delta_pred;
  for (i in 1:N) {
  delta_pred[i][1] = normal_rng(trend[i] * delta_time[i][1], sqrt(psi[i]));
    for (t in 2:T) {
      delta_pred[i][t] = normal_rng(phi[i] * delta[i][t-1] + trend[i] * delta_time[i][t], 
	sqrt(psi[i]));
    }
  }
  
  // obtain log-likelihood (to compute model fit indices)
  vector[N*T] log_lik;
  { // Local environment
    matrix[N,T] temp;
    for (i in 1:N) {
    temp[i][1] = normal_lpdf(delta[i][1] | trend[i] * delta_time[i][1], sqrt(psi[i]));
      for (t in 2:T) {
        temp[i][t] = normal_lpdf(delta[i][t] | phi[i] * delta[i][t-1] + trend[i] * delta_time[i][t],
	 sqrt(psi[i]));
      }
    }
    log_lik = to_vector(temp);
  }
}