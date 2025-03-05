
####################################################################################
#DSEM simulation script (for Chapter 2 of thesis)
#Author: Michael E. Aristodemou
####################################################################################
##############################      Appendix      ##################################
#0. Install & load packages
#1. Generate data from model (DSEM 4 parameters + random effects)
#2. Fit model (DSEM 4 parameters + random effects)
#3. Estimate metrics from Schultzberg & Muthen
#4. Estimate reliability based on Arslan
####################################################################################

####################################################################################
#0. Install & load packages
####################################################################################
start_time <- Sys.time()

#Load packages
library(MASS)
library(dplyr)
library(cmdstanr)
library(beepr)
library(lme4)
library(posterior)
library(bayesplot)
library(ggplot2)
library(kableExtra)
library(tidybayes)
library(beepr)

###################################################################################
#0.1 Top level changes
###################################################################################

# DSEM characteristics
npeople = 1032 # number of subject
tpoints = 20 # number of trials
ndays = 1 # number of days
modtype = "standard" # model type to fit (choice from "standard", "3lvl", "bivariate")

# Number of iterations for simulation
iterations <- 100
####################################################################################
#1. Data generating model (Standard DSEM)
####################################################################################

DataGen <- function(nsubj, nwave, nobs,
                    beta1, beta2, beta3,
                    alpha0, alpha1,
                    sd.subj_tau,
                    sd.wave_loc, sd.wave_scale,
                    cor.subj_tau = matrix(c(1.000, 0.331, 0.223, 0.305, 0.000,
                                            0.331, 1.000, 0.198, 0.282, 0.000,
                                            0.223, 0.198, 1.000, 0.152, 0.000,
                                            0.305, 0.282, 0.152, 1.000, 0.000,
                                            0.000, 0.000, 0.000, 0.000, 1.000),
                                          nrow = 5,
                                          byrow = TRUE)){
  
  ## Generate subject wave and time indicators
  N <- nsubj * nwave * nobs
  subject <- gl(n = nsubj, k = nwave * nobs)
  wave <- rep(gl(n = nwave,k = nobs), nsubj)
  time <- rep(seq_len(nwave*nobs), nsubj)
  mtime <- mean(time)
  
  ## Generate subject, wave level random effects
  # Define correlation matrix
  chol.subj_tau <- t(chol(cor.subj_tau))
  sigma.sd.subj_tau <- diag(sd.subj_tau) %*% chol.subj_tau
  subj_tau <- mvrnorm(nsubj, mu = c(0, 0, 0, 0, 0), 
                      Sigma = sigma.sd.subj_tau %*% t(sigma.sd.subj_tau))
  # Define subject-level random effects
  subj1.tau <- rep(subj_tau[, 1], each = nwave * nobs)
  subj2.tau <- rep(subj_tau[, 2], each = nwave * nobs)
  subj3.tau <- rep(subj_tau[, 3], each = nwave * nobs)
  subj4.tau <- rep(subj_tau[, 4], each = nwave * nobs)
  subj5.tau <- rep(subj_tau[, 5], each = nwave * nobs)
  
  # Define wave level random effects
  wave_loc <- rnorm(nwave, mean = 0, sd = sd.wave_loc)
  wave_scale <- rnorm(nwave, mean = 0, sd = sd.wave_scale)
  wave.loc <- rep(rep(wave_loc, each = nobs), nsubj)
  wave.scale <- rep(rep(wave_scale, each = nobs), nsubj)
  
  # Specify model
  y = rep(0,N)
  y[1] = beta1 
  for(i in 2:N){
    y[i] = rnorm(1, mean = beta1 + subj1.tau[i] + (time[i] - mtime)*(beta2 + subj2.tau[i]) + (y[i-1] - (wave.loc[i] + beta1 + subj1.tau[i]))*(beta3 + subj3.tau[i]) + wave.loc[i], 
                 exp(alpha0 + subj4.tau[i] + (time[i] - mtime)*(alpha1 + subj5.tau[i]) + wave.scale[i]))
  }
  
  # Create data frame
  df.LSME <- data.frame(subject = subject, wave = wave, time = time,
                        y = y)
  
  # Save true values in a vector
  if(modtype == "3lvl") {
    gtruth <- c(beta1, beta2, beta3, 
                cor.subj_tau[1:4,1:4], alpha0,
                sd.wave_loc, sd.subj_tau[1:4])
    }else{
    gtruth <- c(beta1, beta2, beta3, 
              cor.subj_tau[1:4,1:4],  alpha0,
              sd.subj_tau[1:4])
  }
  
  return(list(df.LSME = df.LSME, gtruth = gtruth))
}

#########################################################################
# Simulate data for n repetitions
#########################################################################

fit_results <- list()
generated_data <- list()
ground_truth <- list()

for (i in 1:iterations) {
  
  # Run function
  output <- DataGen(nsubj = npeople, nwave = ndays, nobs = tpoints, # subjects, waves, trials
                    # Fixed effects
                    beta1 = 4.5, # group mean
                    beta2 = 0.05, # group trend 0.000475
                    beta3 = 0.3, # group ar-1
                    alpha0 = 0.5, # group sigma (residual SD (innovations))
                    alpha1 = 0, # group trend on residuals
                    # Random effects
                    sd.subj_tau = c(0.5, 0.01, 0.1, 0.25, 0.0001), # mean sd, trend sd, ar1 sd trend(0.0841)
                    sd.wave_loc = 0.5, # day-to-day sd in means
                    sd.wave_scale = 0, # day-to-daysd in sigma (trialvar)
                    cor.subj_tau = matrix(c(1.000, 0.331, 0.223, 0.305, 0.000,
                                            0.331, 1.000, 0.198, 0.282, 0.000,
                                            0.223, 0.198, 1.000, 0.152, 0.000,
                                            0.305, 0.282, 0.152, 1.000, 0.000,
                                            0.000, 0.000, 0.000, 0.000, 1.000),
                                          nrow = 5,
                                          byrow = TRUE)) # correlation matrix for scale trial-level estimates
  
  
  df.LSME = output$df.LSME
  gtruth = output$gtruth
  # add name labels to gtruth values
  if(modtype == "3lvl"){
  name_vars <- c(
    "alpha", "beta", "phi",
    # Generate rho_u[1,1:10] to rho_u[10,1:10] labels
    with(expand.grid(i = 1:4, j = 1:4), paste0("rho_u[", j, ",", i, "]")),
    "sigma", "tau_d",
    # Generate tau_u[1:10] labels
    paste0("tau_u[", 1:4, "]")
  )
    gtruth[c(2,23)] <- gtruth[c(2,23)] * length(unique(df.LSME$time))
  
  } else if(modtype == "standard"){
    name_vars <- c(
      "alpha", "beta", "phi",
      # Generate rho_u[1,1:10] to rho_u[10,1:10] labels
      with(expand.grid(i = 1:4, j = 1:4), paste0("rho_u[", j, ",", i, "]")),
      "sigma",
      # Generate tau_u[1:10] labels
      paste0("tau_u[", 1:4, "]")
    )
    gtruth[c(2,22)] <- gtruth[c(2,22)] * length(unique(df.LSME$time)) 
  }
  
  
  # create dataframe for gtruth
  gtruth = tibble(
    Parameter = name_vars,
    Value = gtruth
  )
  
  gtruth = gtruth %>% dplyr::rename(".variable" = "Parameter")
  
  # Save generated data and true values for this iteration
  generated_data[[i]] <- df.LSME
  ground_truth[[i]] <- gtruth
  
  ###########################################
  # Preprocessing
  ###########################################
  
  #2-level model
  # Create lagged outcome for AR1
  df.LSME = df.LSME %>% 
    group_by(subject, wave) %>% 
    mutate(Y_tm1 = lag(y)) %>%
    ungroup %>%
    group_by(subject) %>%
    mutate(Y_tm1 = if_else(is.na(Y_tm1), mean(Y_tm1, na.rm=TRUE), Y_tm1)) %>%
    # mutate(Y_tm1 = Y_tm1 - mean(Y_tm1)) %>%
    ungroup()
  
  #3-level model
  df.LSME = df.LSME %>% 
    group_by(subject, wave) %>% 
    mutate(Y_tm1d = lag(y)) %>%
    mutate(Y_tm1d = if_else(is.na(Y_tm1d), mean(Y_tm1d, na.rm=TRUE), Y_tm1d)) %>%
    #mutate(Y_tm1d = Y_tm1d - mean(na.omit(Y_tm1d))) %>%
    ungroup()
  
  # Create within-person centered time
  df.LSME = df.LSME %>%
    group_by(subject) %>%
    mutate(ctime = time/100) %>%
    mutate(mtime = mean(ctime)) %>%
    ungroup()
  
  l2time = df.LSME$mtime[which(df.LSME$time == "1")]  
  
  # Make list of variables and values for DSEM
  dsem_list <- list(N_subj = length(unique(as.numeric(as.factor(df.LSME$subject)))), # subject number
                    Y = df.LSME$y, # outcome variable matrix
                    N_obs = length(unique(df.LSME$time)),
                    N_days = length(unique(as.numeric(as.factor(df.LSME$wave)))),
                    subj = as.numeric(as.factor(df.LSME$subject)),
                    day = as.numeric(as.factor(df.LSME$wave)),
                    N = nrow(df.LSME),
                    mtime = l2time,
                    time = df.LSME$ctime,
                    Y_tm1 = df.LSME$Y_tm1, # wp centered AR1
                    Y_tm1d = df.LSME$Y_tm1d) # wp/wd centered AR1
  
  # View data structure
  #str(dsem_list) 
  
  # Quick visualization of data for 1 person
  #plot(df.LSME$y[1:(length(unique(df.LSME$time)))],type="l") + # lineplot
  #  abline(h=mean(df.LSME$y[1:(length(unique(df.LSME$time)))], na.rm=TRUE),col="dodgerblue2") # mean line
  
  #########################################################################################################
  #2. Fit model to data
  #########################################################################################################
  
  #Compile model
  stan_spec = paste0("~\\Simulations_DSEM\\dsem_", modtype, ".stan")
  model2 <- cmdstan_model(stan_file = stan_spec, stanc_options = list("O1"))
  # Sample from the model
  # dsem_list$grainsize <- 1 #set grainsize (how many subjects per thread)
  time2 = system.time(fit2 <- model2$sample(dsem_list,
                                            chains = 4,
                                            parallel_chains = 4,
                                            iter_warmup = 1000,
                                            iter_sampling = 2000,
                                            save_warmup = FALSE))
  
  # Save the fit result for this iteration
  fit_results[[i]] <- fit2
  
  # Optionally print progress
  cat("Iteration", i, "completed\n")
}

# Save results if needed
# Create the filename dynamically using paste0
fit_name <- paste0("~\\Simulations_DSEM\\model_fit\\fit_results_n", npeople, "t", tpoints, "w", ndays, ".RData")
sim_name <- paste0("~\\Simulations_DSEM\\sim_data\\sim_data_n", npeople, "t", tpoints, "w", ndays, ".RData")
save(fit_results, file = fit_name) # use paste0 to automatically set iter
save(generated_data, ground_truth, file = sim_name)

#######################################################################################################
#3. Extract indices from model fit object
#######################################################################################################

# Define the folder path
folder_path <- "~/Simulations_DSEM/model_fit/fits/"
# Delete all files in the folder
unlink(paste0(folder_path, "*"))

##############################################################
#3.1 Create data frame with all variables to estimate metrics
##############################################################

# Define file paths to save intermediate results
save_file <- "~\\Simulations_DSEM\\model_fit\\fits\\est_list_progress_iter_"
progress_file <- "~\\Simulations_DSEM\\model_fit\\fits\\iteration_progress.txt"

# Initialize or load progress
if (file.exists(progress_file)) {
  start_iteration <- as.numeric(readLines(progress_file))
} else {
  start_iteration <- 1  # Start from scratch
}

# Start the loop from where it stopped
for (i in start_iteration:iterations) {
  # Try-catch block to handle out of memory or other issues
  tryCatch({
    # Estimate results for the current iteration
    current_result <- fit_results[[i]] %>%
      gather_draws(`alpha`, `beta`, `phi`,
                   `sigma`, `tau_.+`, `rho_.+`,
                   regex = TRUE) %>%
      ungroup() %>%
      group_by(.variable) %>%
      summarise(mean = mean(.value),
                se = sd(.value),
                ci_lower = quantile(.value, 0.025),
                ci_upper = quantile(.value, 0.975)) %>%
      ungroup()
    
    # Bind current results to gtruth
    current_result = left_join(current_result, gtruth, by = ".variable") %>%
      rename("gtruth" = "Value")
    
    # Save the current result immediately after calculation
    saveRDS(current_result, paste0(save_file, i, ".rds"))
    
    # Save progress after each iteration
    writeLines(as.character(i), progress_file)
    
    # Clear the current result from memory
    rm(current_result)
    
    # Perform garbage collection every 3 iterations
    if (i %% 1 == 0) {
      message("Running garbage collection (gc()) after iteration: ", i)
      gc()
    }
    
  }, error = function(e) {
    # Save current iteration in case of error
    writeLines(as.character(i), progress_file)
    
    # Print error message and stop the loop
    message("Error at iteration: ", i)
    message("Error message: ", e$message)
    stop("Stopping loop due to memory issue.")
  })
}

# Load all saved results
est_list <- list()

for (i in 1:iterations) {
  result_file <- paste0("~\\Simulations_DSEM\\model_fit\\fits\\est_list_progress_iter_", i, ".rds")
  if (file.exists(result_file)) {
    est_list[[i]] <- readRDS(result_file)
  }
}


# Estimate difference between mean estimate and true score
for (i in 1:iterations) {
  est_list[[i]] = est_list[[i]] %>%
    group_by(.variable) %>%
    mutate(abs_bias = (mean - gtruth)^2) %>%
    mutate(inc95CI = if_else(gtruth >= ci_lower & gtruth <= ci_upper, 1, 0)) %>%
    mutate(nonzeroCI = if_else(ci_lower <= 0 & ci_upper >= 0, 0, 1)) %>%
    ungroup()
}

# Save the list to disk
list_name <- paste0("~\\Simulations_DSEM\\model_fit\\est_list_", npeople, "t", tpoints, "w", ndays, ".rds")
saveRDS(est_list, file = list_name)
# est_list = readRDS(list_name)


# Estimate mean and sd over sim repetitions
est_df <- bind_rows(est_list)  # Converts list to a single data frame

# Create a df with all data needed for metrics
rep_list = est_df %>%
  group_by(.variable) %>%
  summarise(mean_est = mean(mean),
            mean_se = mean(se),
            sd_est = sd(mean),
            sum_bias = sum(abs_bias),
            sum_incCI = sum(inc95CI),
            sum_nozeroCI = sum(nonzeroCI)
  ) %>%
  ungroup()

rep_list = left_join(rep_list, gtruth, by = ".variable") %>%
  rename("gtruth" = "Value")

##################################################
# Create function to estimate metrics
##################################################

# this function estimates all the metrics
metrics = function(data){
  metrics = data %>%
    group_by(.variable) %>%
    mutate(rel_bias = mean_est/gtruth) %>%
    mutate(se_sd = mean_se/sd_est) %>%
    mutate(mse = sum_bias/iterations) %>%
    mutate(cov95 = sum_incCI/iterations) %>%
    mutate(power = sum_nozeroCI/iterations) %>%
    ungroup()
  
  return(metrics)
}

out = metrics(data = rep_list)
# save to csv
power_name = paste0("~\\Simulations_DSEM\\power\\power_", modtype, "_n", npeople, "t", tpoints, "w", ndays, ".csv")
write.csv(out, file = power_name, row.names = FALSE)

rm(est_list, est_df, rep_list)
gc()

####################################################
# Reliability function
####################################################

dsem_reliability <- function(fit_model) {
  
  # Get the SD of estimates (true score variance)
  true_score_var1 <- fit_model %>%
    gather_draws(`tau_.+`, regex = TRUE) %>% 
    mutate(.variable = stringr::str_sub(.variable, 1)) %>% 
    mutate(latent_var = .value^2) %>% 
    group_by(.variable) %>%
    summarise(latent_var = mean(latent_var))
  
  # Get standard error of estimates (SEE2)
  if(modtype == "3lvl"){
    SEE2_1 <- fit_model %>%
      gather_draws(`u`[id, parameter],`d`[id], regex = TRUE) %>%
      ungroup() %>%
      tidyr::unite(.variable, .variable, parameter, sep = "__") %>%
      group_by(.variable, id) %>%
      summarise(SEE2 = var(.value)) %>%
      summarise(SEE2 = mean(SEE2))
  } else {
    SEE2_1 <- fit_model %>%
      gather_draws(`u`[id, parameter], regex = TRUE) %>%
      ungroup() %>%
      tidyr::unite(.variable, .variable, parameter, sep = "__") %>%
      group_by(.variable, id) %>%
      summarise(SEE2 = var(.value)) %>%
      summarise(SEE2 = mean(SEE2))
  }
  
  # Create a table with true score variance and SEE2
  rel_table <- true_score_var1
  rel_table$SEE2 <- SEE2_1$SEE2
  
  # Estimate reliability
  rel_table <- rel_table %>%
    mutate(reliability = 1 - SEE2 / latent_var)
  
  # Rename the .variable column (add ifelse here)
  if(modtype == "3lvl"){
    rel_table <- rel_table %>%
      mutate(.variable = case_when(
        .variable == "tau_u[1]" ~ "alpha",
        .variable == "tau_u[2]" ~ "beta",
        .variable == "tau_u[3]" ~ "phi",
        .variable == "tau_u[4]" ~ "sigma",
        .variable == "tau_d" ~ "day_sd",
        TRUE ~ .variable
    ))
  }
  else{
    rel_table <- rel_table %>%
      mutate(.variable = case_when(
        .variable == "tau_u[1]" ~ "alpha",
        .variable == "tau_u[2]" ~ "beta",
        .variable == "tau_u[3]" ~ "phi",
        .variable == "tau_u[4]" ~ "sigma",
        TRUE ~ .variable
      ))
  }
  
  return(rel_table)
}

####################################################
# Run reliability function
####################################################

# Define file paths to save intermediate results
save_file <- "~\\Simulations_DSEM\\model_fit\\fits\\reliability_list_progress_iter_"
progress_file <- "~\\Simulations_DSEM\\model_fit\\fits\\reliability_iteration_progress.txt"

# Initialize or load progress
if (file.exists(progress_file)) {
  start_iteration <- as.numeric(readLines(progress_file))  # Resume from last saved iteration
} else {
  start_iteration <- 1  # Start from the beginning
}

# Start the loop from where it stopped
for (i in start_iteration:iterations) {
  # Try-catch block to handle out-of-memory or other issues
  tryCatch({
    # Compute reliability for the current iteration
    current_reliability <- dsem_reliability(fit_model = fit_results[[i]])
    
    # Save the current result immediately after calculation
    saveRDS(current_reliability, paste0(save_file, i, ".rds"))
    
    # Save progress after each iteration
    writeLines(as.character(i), progress_file)
    
    # Clear the current result from memory
    rm(current_reliability)
    
    # Perform garbage collection every 3 iterations
    if (i %% 3 == 0) {
      message("Running garbage collection (gc()) after iteration: ", i)
      gc()
    }
    
  }, error = function(e) {
    # Save progress in case of an error
    writeLines(as.character(i), progress_file)
    
    # Print error message and stop the loop
    message("Error at iteration: ", i)
    message("Error message: ", e$message)
    stop("Stopping loop due to memory issue.")
  })
}

# Load all saved reliability results
reliability_list <- list()

for (i in 1:iterations) {
  result_file <- paste0("~\\Simulations_DSEM\\model_fit\\fits\\reliability_list_progress_iter_", i, ".rds")
  if (file.exists(result_file)) {
    reliability_list[[i]] <- readRDS(result_file)
  }
}

reliability_df <- bind_rows(reliability_list)  # Converts list to a single data frame

# Create a df with all data needed for metrics
rel_df = reliability_df %>%
  group_by(.variable) %>%
  summarise(mean_rel = mean(reliability),
            sd_rel = sd(reliability),
            lower_rel = quantile(reliability, 0.025),
            upper_rel = quantile(reliability, 0.975),
            # SEE2
            mean_see2 = mean(SEE2),
            sd_see2 = sd(SEE2),
            lower_see2 = quantile(SEE2, 0.025),
            upper_see2 = quantile(SEE2, 0.975),
            # latent var
            mean_lvar = mean(latent_var),
            sd_lvar = sd(latent_var),
            lower_lvar = quantile(latent_var, 0.025),
            upper_lvar = quantile(latent_var, 0.975),
  ) %>%
  ungroup()

# save reliability metrics to csv
rel_name = paste0("~\\Simulations_DSEM\\reliability\\reliability_", modtype, "_n", npeople, "t", tpoints, "w", ndays, ".csv")
write.csv(rel_df, file = rel_name, row.names = FALSE)

#########################################
# log time it takes for script
#########################################

# Log the end time
end_time <- Sys.time()

# Calculate the duration and print it
total_time_taken <- end_time - start_time
print(paste("Total execution time:", total_time_taken))

