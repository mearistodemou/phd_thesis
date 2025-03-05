###################
# Chapter 2 tables
###################

library(knitr)
library(kableExtra)
library(tidyverse)
library(tidyr)

################################
# Load dataset
################################

table = read.csv("~\\Simulations_DSEM\\simulation_estimates.csv")

# Preprocess data
result <- table %>%
  select(gtruth, .variable, model) %>%
  mutate(row_id = row_number()) %>%  # Add a row identifier
  pivot_wider(names_from = model, 
              values_from = gtruth, 
              names_prefix = "model") %>%  # Spread the model column
  select(-row_id)  # Remove the row identifier

# Drop rows with NA in each model column
mod1 <- result %>%
  select(.variable, model1) %>%
  drop_na()

mod2 <- result %>%
  select(.variable, model2) %>%
  drop_na()

mod3 <- result %>%
  select(.variable, model3) %>%
  drop_na()

mod4 <- result %>%
  select(.variable, model4) %>%
  drop_na()

# Join the four dataframes by `.variable`
table <- mod1 %>%
  full_join(mod2, by = ".variable") %>%
  full_join(mod3, by = ".variable") %>%
  full_join(mod4, by = ".variable")

# Order by variable names
row_order = c("alpha1_mean", "beta1_mean", "phi1_mean", "sigma1_mean",
              "gamma1_mean", "gamma2_mean",
              "alpha2_mean", "beta2_mean", "phi2_mean", "sigma2_mean",
              # scale model
              "alpha1_sd", "beta1_sd", "phi1_sd", "sigma1_sd",
              "gamma1_sd", "gamma2_sd",
              "alpha2_sd", "beta2_sd", "phi2_sd", "sigma2_sd",
              "day_sd")


table <- table %>%
  mutate(.variable = factor(.variable, levels = row_order)) %>%
  arrange(.variable)

#################################
# Table for simualtion estimates
#################################

#Output latex code apa table
#Unstandardized parameters
kable(
  table, #subset to desired columns
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Ground truth for model parameters", #change caption
  col.names = c("", "Model1", "Model 2", "Model 3", "Model 4"), #rename your columns
  align = "c", #center columns
) %>%
  row_spec(row = 0, align = "c")%>%
  add_footnote(c("Note: PHI = inertia, alpha = mean",
                 "sigma = trial-to-trial variability variability, beta = trend",
                 "gamma = coupling", "day_sd = day-to-day variability"))
