############################################
# DSEM for Chapter 6 of PhD thesis
# Methylphenidate modulates RTV
############################################
# Author: Michael E. Aristodemou
#############   Appendix    ################
# 0. Load packages
# 1. Load dataset
# 2. Preprocessing
# 3. Sample from DSEM
# 4. Model diagnostics
############################################

############################################
# 0. Load packages
############################################

{
  #library(cmdstanr)
  library(dplyr)
  library(bayesplot)
  library(tidyr)
  library(loo)
  library(rstantools)
  library(ggplot2)
  library(rstanarm)
  library(beepr)
  library(corrplot)
  library(rstan)
}

############################################
# 1. Load dataset
############################################

# Choose drug condition
drug_con <- "MPH" # Options: "MPH", "SUL", "PBO"

# Load datafile
data_string <- paste0("~\\CatVar\\Data\\df_", drug_con, ".csv")
data <- read.csv(data_string)

# Select model
mod_ext = "~\\CatVar\\Chapter6\\VarEstimates\\DSEM\\Stan\\Models\\"
dsem_stand = paste0(mod_ext, "dsem_standard.stan")
dsem_surp = paste0(mod_ext, "dsem_surprise.stan")
############################################
# 2. Preprocessing
############################################

# Count the number of trials for each subject
trial_counts <- data %>%
  group_by(subject) %>%
  summarise(t_count = n())

# Get SubjectIDs with 56 or more trials
valid_subjects <- trial_counts %>%
  filter(t_count >= 56) %>%
  pull(subject)

# Filter the original data to include only valid subjects
data <- data %>%
  filter(subject %in% valid_subjects)

# Replace missing logRT values with a marker (999)
na_indices <- is.na(data$logRT)
data$logRT[na_indices] <- 999

# Function to count the missing markers (999)
computeMisses <- function(x) {
  sum(x == 999)
}

# Determine unique subjects and number of trials per subject
unique_subjects <- unique(data$subject)
N_sub <- length(unique_subjects)
T_trials <- length(unique(data$trial))  # assumes equal number of trials per subject

# Create matrices for outcome (Y) and time
Y_mat <- matrix(data$logRT, nrow = N_sub, byrow = TRUE)
time_mat <- matrix(data$trial, nrow = N_sub, byrow = TRUE)
sup_mat <- matrix(data$surprise, nrow = N_sub, byrow = TRUE)

# Compute coordinates for imputation (for every occurrence of 999 in Y_mat)
coordinates <- NULL
for (it in 1:nrow(Y_mat)) {
  miss_idx <- which(Y_mat[it, ] == 999)
  if (length(miss_idx) > 0) {
    coordinates <- rbind(coordinates, cbind(it, miss_idx))
  }
}

# Build the DSEM list
dsem_list <- list(
  subid  = unique_subjects,
  N      = N_sub,
  Y      = Y_mat,
  time   = time_mat,
  surp   = sup_mat,
  T      = T_trials,
  N_miss = computeMisses(data$logRT),
  ii_miss = coordinates,
  y_mean = mean(data$logRT[data$logRT != 999]),
  y_sd   = sd(data$logRT[data$logRT != 999])
)

str(dsem_list)


#################################################
# 3. Sample from DSEM
#################################################

# rstan sampler
mod <- stan(file = dsem_stand, # choose your model
            data = dsem_list, # insert data
            verbose = FALSE, # limited sampling output
            iter = 2000, chains = 4, # MCMC details
            cores = 4) # PC cores (1 for each chain)

# summarize estimates
sum = summary(mod)

#################################################
# 4. Model diagnostics
#################################################

# Save posterior draws
mod_draw = rstan::extract(mod)
mod_draw = as.array(mod)
np_cp <- nuts_params(mod)

# 5.1 Plot rhats
rhats = rhat(mod)
rhats_sub = rhats[grepl("gamma|tau_Omega", names(rhats))]

hatplot = mcmc_rhat(rhats_sub) + yaxis_text(hjust = 1)

# 5.2 Plots traceplots
traceplot = mcmc_trace(mod_draw, pars = c("gamma[1]", "gamma[2]", "gamma[3]", "gamma[4]",
                              "tau_Omega[1]", "tau_Omega[2]", "tau_Omega[3]", 
                              "tau_Omega[4]"), np = np_cp) + xlab("Post-warmup iteration")

# 5.3 Autocorrelation
acfplot = mcmc_acf(mod_draw, pars = c("gamma[1]", "gamma[2]", "gamma[3]", "gamma[4]",
                                       "tau_Omega[1]", "tau_Omega[2]", "tau_Omega[3]", 
                                       "tau_Omega[4]"), lags = 10)

# 5.4 Histograms
areaplot = mcmc_areas(mod_draw,
           pars = c("gamma[1]", "gamma[2]", "gamma[3]", "gamma[4]",
                    "tau_Omega[1]", "tau_Omega[2]", "tau_Omega[3]", 
                    "tau_Omega[4]"),
           prob = 0.80)

# 5.5 Name plot files
hat_name <- paste0("~\\CatVar\\Stan\\Visualization\\rhat_", drug_con, ".png")
trace_name <- paste0("~\\CatVar\\Stan\\Visualization\\trace_", drug_con, ".png")
hist_name <- paste0("~\\CatVar\\Stan\\Visualization\\hist_", drug_con, ".png")
acf_name <- paste0("~\\CatVar\\Stan\\Visualization\\acf_", drug_con, ".png")

# 5.6 Save plots
ggsave(plot = hatplot, hat_name, dpi = 300, height = 10, width = 15, bg = "white")
ggsave(plot = traceplot, trace_name, dpi=300, height = 10, width = 15, bg = "white")
ggsave(plot = areaplot, hist_name, dpi=300, height = 10, width = 15, bg = "white")
ggsave(plot = acfplot, acf_name, dpi=300, height = 10, width = 15, bg = "white")
