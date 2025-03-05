#################################
# Chapter 7: Stan output tables
# Author: Michael E. Aristodemou
#################################

library(knitr)
library(kableExtra)
library(tidyverse)
library(dplyr)

#######################################
#Select data
#######################################

table = as.data.frame(sum1) # for pc model
table = as.data.frame(sum2) # for vp model
table = as.data.frame(sum3) # for ss model

# Remove random effect covariances (these are plotted)
table <- table %>%
  filter(!grepl("^rho_u", variable))


# Use kable without specifying row.names
kable(
  table, # Subset to desired columns if needed
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Variable Precision model estimates", # Change caption
  col.names = c("Parameter", "Mean Est.", "Median Est.",
                "SD Est.", "MAD", "Q5",
                "Q95", "Rhat", "ESS bulk", "ESS tail"), # Rename your columns
  align = "c" # Center columns
) %>%
  row_spec(row = 0, align = "c")



