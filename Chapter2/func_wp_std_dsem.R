#####################################################
# Estimate within-person standardized coefficients
# Using method from Schuurman et al., 2016 (Mplus)

# Michael E. Aristodemou
#####################################################

################
# load packages
################
library(tidyr)
library(dplyr)

######################################################################
# Create function to calculate within-person standardized parameters
######################################################################

wp_std_effects <- function(dataset, fit, outcome, lag1, time, sub_id) {
  # Estimate SD of predictors and outcome
  df_sub = dataset %>%
    group_by(!! sym(sub_id)) %>%
    summarise(
      ySD = sd(!! sym(outcome)),
      lag1SD = sd(!! sym(lag1)),
      timeSD = sd(!! sym(time))
    ) %>%
    ungroup()
  
  # Extract fixed-effect parameter estimates (median of posterior)
  b1_draws = fit$draws(format = "df", variables = c("beta", "phi")) %>%
    summarise(
      beta = median(beta),
      phi = median(phi)
    )
  
  # Extract subject-specific random effects
  subj_effects <- fit %>%
    gather_draws(`u`[id, parameter], regex = TRUE) %>%
    ungroup() %>%
    group_by(id, parameter) %>%
    summarise(u = median(.value)) %>%
    ungroup() %>%
    pivot_wider(
      names_from = parameter,
      values_from = u,
      names_prefix = "u_"
    ) %>%
    mutate(id = as.factor(id))
  
  # Merge subject-specific effects with subject-specific SDs
  wp_effects = left_join(df_sub, subj_effects, by = c("subject" = "id")) %>%
    mutate(
      beta_wp = b1_draws$beta + u_2,
      phi_wp = b1_draws$phi + u_3,
      phi_stdw = phi_wp * (lag1SD / ySD),
      beta_stdw = beta_wp * (timeSD / ySD)
    )
  
  # Calculate mean within-person standardized effects
  results <- list(
    phi_stdw_mean = mean(wp_effects$phi_stdw, na.rm = TRUE),
    beta_stdw_mean = mean(wp_effects$beta_stdw, na.rm = TRUE),
    wp_effects = wp_effects
  )
  
  return(results)
}

###################
# Run function
###################

wp_std_effects(dataset = df.LSME, # add the dataset containing the outcome, the lagged predictor, and the trend predictor
               fit = fit2, # add the Stan fit object, e.g., fit2
               outcome = "y", # add the name of the outcome column, e.g., "y"
               lag1 = "Y_tm1", # add the name of the lagged predictor column, e.g., "Y_tm1"
               time = "ctime", # add the name of the time predictor column, e.g., "ctime"
               sub_id = "subject" # add the name of the subject id column, e.g., "subject"
               ) # you get a list called results
