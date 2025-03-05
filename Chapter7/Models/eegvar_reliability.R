#################################
# EEGVar Reliability
#################################

library(tidybayes)
library(brms)
library(cmdstanr)
library(MASS)
library(dplyr)
library(beepr)
library(lme4)
library(posterior)
library(bayesplot)
library(ggplot2)
library(kableExtra)

####################################################
# Choose regional vs global model
####################################################

regions = "global"
regions = "posterior"

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
  
    SEE2_1 <- fit_model %>%
      gather_draws(`u`[id, parameter], regex = TRUE) %>%
      ungroup() %>%
      tidyr::unite(.variable, .variable, parameter, sep = "__") %>%
      group_by(.variable, id) %>%
      summarise(SEE2 = var(.value)) %>%
      summarise(SEE2 = mean(SEE2))
  
  # Create a table with true score variance and SEE2
  rel_table <- true_score_var1
  rel_table$SEE2 <- SEE2_1$SEE2
  
  # Estimate reliability
  rel_table <- rel_table %>%
    mutate(reliability = 1 - SEE2 / latent_var)
  
  # Rename the .variable column (add ifelse here)
    rel_table <- rel_table %>%
      mutate(.variable = case_when(
        .variable == "tau_u[1]" ~ "alpha1",
        .variable == "tau_u[2]" ~ "beta1",
        .variable == "tau_u[3]" ~ "phi1",
        .variable == "tau_u[4]" ~ "couple1",
        .variable == "tau_u[5]" ~ "sigma1",
        .variable == "tau_u[6]" ~ "alpha2",
        .variable == "tau_u[7]" ~ "beta2",
        .variable == "tau_u[8]" ~ "phi2",
        .variable == "tau_u[9]" ~ "couple3",
        .variable == "tau_u[10]" ~ "sigma2",
        
        TRUE ~ .variable
      ))

  
  return(rel_table)
}

####################################################
# Run reliability function
####################################################

# Compute reliability for the current iteration
rel_df <- dsem_reliability(fit_model = fit3)
    
# save reliability metrics to csv
rel_name = paste0("~\\EEGvar\\Visualization\\Diagnostics\\reliability_ss_", regions, ".csv")
write.csv(rel_df, file = rel_name, row.names = FALSE)
