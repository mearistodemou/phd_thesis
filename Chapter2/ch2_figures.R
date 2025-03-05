###############################################
# Figures for Chapter 2: Analytical framework
###############################################

#############################################
# Visualize simulations results
# Author: Michael E. Aristodemou
# email: michael.aristodemou@radboudumc.nl
#############################################

#####################################
########   Appendix    ##############
#0. Load packages & make decisions
#1. Load datasets
#2. Plot posteriors
#3. Plot relative bias
#4. Plot SE/SD
#5. Plot MSE
#6. Plot 95% coverage
#7. Plot power
#8. Plot reliability
######################################

#0. load packages
library(ggplot2)
library(dplyr)
library(viridis)

#############################
# Top-level changes
#############################

###################################################
# Choose one of four options
# Or enter your own values + modchoice = "custom
# Then just ctrl+A --> ctrl+Enter
# To get all figures
###################################################


#####################
#1. Load data
#####################

# retrive filename

power_name_std1 <- paste0("~\\Simulations_DSEM\\power\\power_standard_n100t50w1.csv")
power_name_std2 <- paste0("~\\Simulations_DSEM\\power\\power_standard_n1032t20w1.csv")
power_name_3lvl <- paste0("~\\Simulations_DSEM\\power\\power_3lvl_n100t10w20.csv")
power_name_bivr <- paste0( "~\\Simulations_DSEM\\power\\power_bivariate_n150t100w1.csv")

#1. load datasets
est1 = read.csv(power_name_std1)
est2 = read.csv(power_name_std2)
est3 = read.csv(power_name_3lvl)
est4 = read.csv(power_name_bivr)

# Model id
est1$model <- 1
est2$model <- 2
est3$model <- 3
est4$model <- 4

# Bind all dataframes
est_df = rbind(est1,est2,est3,est4)

# Reorder values
var_order = c("alpha1_mean", "beta1_mean", "phi1_mean", "sigma1_mean", "gamma1_mean",
              "alpha2_mean", "beta2_mean", "phi2_mean", "sigma2_mean", "gamma2_mean",
              "alpha1_sd", "beta1_sd", "phi1_sd", "sigma1_sd", "gamma1_sd",
              "alpha2_sd", "beta2_sd", "phi2_sd", "sigma2_sd", "gamma2_sd",
              "day_sd",
              with(expand.grid(i = 1:10, j = 1:10), paste0("rho_u[", j, ",", i, "]"))
)

# rename variables
est_df <- est_df %>%
  mutate(.variable = ifelse(.variable == "alpha", "alpha1", 
                            ifelse(.variable == "phi", "phi1",
                                   ifelse(.variable == "beta", "beta1",
                                          ifelse(.variable == "sigma", "sigma1", .variable)))))
# rename final
est_df <- est_df %>%
  mutate(.variable = case_when(
    .variable == "alpha1" ~ "alpha1_mean", 
    .variable == "beta1" ~ "beta1_mean",
    .variable == "phi1" ~ "phi1_mean",
    .variable == "sigma1" ~ "sigma1_mean",
    .variable == "couple1" ~ "gamma1_mean",
    .variable == "alpha2" ~ "alpha2_mean", 
    .variable == "beta2" ~ "beta2_mean",
    .variable == "phi2" ~ "phi2_mean",
    .variable == "sigma2" ~ "sigma2_mean",
    .variable == "couple3" ~ "gamma2_mean",
    # random effects
    .variable == "tau_u[1]" ~ "alpha1_sd", 
    .variable == "tau_u[2]" ~ "beta1_sd",
    .variable == "tau_u[3]" ~ "phi1_sd",
    .variable == "tau_u[4]" ~ "sigma1_sd",
    .variable == "tau_u[5]" ~ "gamma1_sd",
    .variable == "tau_u[6]" ~ "alpha2_sd",
    .variable == "tau_u[7]" ~ "beta2_sd",
    .variable == "tau_u[8]" ~ "phi2_sd",
    .variable == "tau_u[9]" ~ "sigma2_sd",
    .variable == "tau_u[10]" ~ "gamma2_sd",
    .variable == "tau_d" ~ "day_sd",
    TRUE ~ .variable  # Keep other values unchanged
  )) %>%
  mutate(
    .variable = factor(.variable, levels = rev(var_order)))
  

# Remove rows where .variable starts with "rho_"
est_df <- est_df %>%
  filter(!grepl("^rho_", .variable))

write.csv(est_df, "~\\Simulations_DSEM\\simulation_estimates.csv", row.names = FALSE)

###########################################################################################
# 2. Plot relative bias
###########################################################################################

rel_tile <- ggplot(est_df, aes(x = model, y = .variable, fill = rel_bias)) + 
  geom_tile(color = "black", size = 1.5) +
  geom_text(aes(label = round(rel_bias, 2)), color = "white", size = 4) +  # Add values
  scale_fill_viridis_c(
    name = "Relative Bias",  # Add a legend title
    limits = c(0, 2),        # Set the range for the fill values
    option = "viridis"       # Default Viridis color palette
  ) +
  theme_minimal() +  # Optional: Clean minimalistic theme
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )+
  labs(x = "Model Type")# Optional: Clean minimalistic theme

# Color coding y-axis labels based on specific .variable values

# Save figure
ggsave(plot = rel_tile,
       file = "~\\Simulations_DSEM\\visualizations\\dsem_tile_relbias.png",
       dpi = 300,
       width = 6,
       height = 10,
       bg = "white")


###########################################################################################
# 3. Plot SE/SD
###########################################################################################

se_sd_tile <- ggplot(est_df, aes(x = model, y = .variable, fill = se_sd)) + 
  geom_tile(color = "black", size = 1.5) +
  geom_text(aes(label = round(se_sd, 2)), color = "white", size = 4) +  # Add values
  scale_fill_viridis_c(
    name = "SE/SD",  # Add a legend title
    #limits = c(0, 2),        # Set the range for the fill values
    oob = scales::squish,    # Squish values outside the range
    option = "viridis"       # Default Viridis color palette
  ) +
  theme_minimal() +  # Optional: Clean minimalistic theme
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )+
  labs(x = "Model Type")# Optional: Clean minimalistic theme

# Color coding y-axis labels based on specific .variable values

# Save figure
ggsave(plot = se_sd_tile,
       file = "~\\Simulations_DSEM\\visualizations\\dsem_tile_se_sd.png",
       dpi = 300,
       width = 6,
       height = 10,
       bg = "white")

###########################################################################################
# 4. MSE
###########################################################################################

mse_tile <- ggplot(est_df, aes(x = model, y = .variable, fill = mse)) + 
  geom_tile(color = "black", size = 1.5) +
  geom_text(aes(label = round(mse, 2)), color = "white", size = 4) +  # Add values
  scale_fill_viridis_c(
    name = "MSE",  # Add a legend title
    #limits = c(0, 2),        # Set the range for the fill values
    #oob = scales::squish,    # Squish values outside the range
    option = "viridis"       # Default Viridis color palette
  ) +
  theme_minimal() +  # Optional: Clean minimalistic theme
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )+
  labs(x = "Model Type")# Optional: Clean minimalistic theme

# Color coding y-axis labels based on specific .variable values

# Save figure
ggsave(plot = mse_tile,
       file = "~\\Simulations_DSEM\\visualizations\\dsem_tile_mse.png",
       dpi = 300,
       width = 6,
       height = 10,
       bg = "white")

###########################################################################################
# 5. 95% CI
###########################################################################################

cov95_tile <- ggplot(est_df, aes(x = model, y = .variable, fill = cov95)) + 
  geom_tile(color = "black", size = 1.5) +
  geom_text(aes(label = round(cov95, 2)), color = "white", size = 4) +  # Add values
  scale_fill_viridis_c(
    name = "95% Coverage",  # Add a legend title
    #limits = c(0, 2),        # Set the range for the fill values
    oob = scales::squish,    # Squish values outside the range
    option = "viridis"       # Default Viridis color palette
  ) +
  theme_minimal() +  # Optional: Clean minimalistic theme
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )+
  labs(x = "Model Type")# Optional: Clean minimalistic theme

# Color coding y-axis labels based on specific .variable values

# Save figure
ggsave(plot = cov95_tile,
       file = "~\\Simulations_DSEM\\visualizations\\dsem_tile_cov95.png",
       dpi = 300,
       width = 6,
       height = 10,
       bg = "white")

###########################################################################################
# 6. Power
###########################################################################################

power_tile <- ggplot(est_df, aes(x = model, y = .variable, fill = power)) + 
  geom_tile(color = "black", size = 1.5) +
  geom_text(aes(label = round(power, 2)), color = "white", size = 4) +  # Add values
  scale_fill_viridis_c(
    name = "Power",  # Add a legend title
    #limits = c(0, 2),        # Set the range for the fill values
    oob = scales::squish,    # Squish values outside the range
    option = "viridis"       # Default Viridis color palette
  ) +
  theme_minimal() +  # Optional: Clean minimalistic theme
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )+
  labs(x = "Model Type")# Optional: Clean minimalistic theme

# Color coding y-axis labels based on specific .variable values

# Save figure
ggsave(plot = power_tile,
       file = "~\\Simulations_DSEM\\visualizations\\dsem_tile_power.png",
       dpi = 300,
       width = 6,
       height = 10,
       bg = "white")



# Calculate the median of se_sd when model == 4
sum_stat <- est_df %>%
  summarise(sum_stat = median(cov95, na.rm = TRUE)) %>%
  pull(sum_stat)

# Print the result
print(sum_stat)


# Calculate the median of se_sd when model == 4
sum_stat <- est_df %>%
  filter(model == 1) %>%
  summarise(sum_stat = median(mse, na.rm = TRUE)) %>%
  pull(sum_stat)

# Print the result
print(sum_stat)




