##########################################################
# A sequence sensitive model of encoding precision
# Model specification, estimation, and comparison script

# Author: Michael E. Aristodemou
##########################################################

#######################################
##########    Appendix    #############
#0. Load data and packages
#1. Population coding model
#2. Variable precision model
#3. Sequence sensitive model
#4. Compare all models using loo
#######################################

#######################################
# Top-level changes
#######################################

#regions = "posterior"
regions = "global"

#######################################
#0. Load data and packages
#######################################

#0.1 Load packages
{
library(cmdstanr)
library(dplyr)
library(bayesplot)
library(tidyr)
library(loo)
library(rstantools)
library(ggplot2)
library(rstanarm)
library(beepr)
library(corrplot)
}

#0.2 Load data
{
if (regions == "global"){
eegdat <- read.csv("~/EEGvar/Data/eegnoise_allpp.csv") # whole brain
} else if (regions == "posterior"){
eegdat <- read.csv("~/EEGvar/Data/eegnoise_posterior.csv") # posterior parietal & occipital regions
}

RTdat <- read.csv("~/EEGvar/Behaviour/eegRT.csv") # response times
}

#0.3 Restructuring data
# Create item counter
eegdat<-eegdat %>%
  dplyr::group_by(subject) %>%
  dplyr::arrange(subject,epoch) %>%
  dplyr::mutate(grp = cumsum(row_number(1))) %>%
  dplyr::group_by(subject, trialnr, grp) %>%
  dplyr::mutate(item_count = row_number()) %>%
  ungroup() %>%
  select(-grp)

# Create trial level 1/f int and slope
# Offset (intercept)
eegdat <- eegdat %>%
  group_by(subject, trialnr) %>%
  mutate(trial_offset = mean(item_offset, na.rm = TRUE)) %>%
  ungroup
# Exponent (slope)
eegdat <- eegdat %>%
  group_by(subject, trialnr) %>%
  mutate(trial_exponent = mean(item_pink, na.rm = TRUE)) %>%
  ungroup

# Create trial level eegdat
# Remove all rows with an item counter > 1 (remove duplicate trials)
eegtrial = eegdat[which(eegdat$item_count == "1"),]

# Load behavioral data
jointdat<-left_join(eegtrial, RTdat, join_by(subject==Subject, trialnr==TrialNum))
# remove subject 4 (this is a temporary fix)
jointdat <- jointdat[!(jointdat$subject %in% c("4","34")),]
#Create trial_noise estimate to predict RT

# Substitute NA with person mean
jointdat <- jointdat %>% 
  group_by(subject) %>% 
  mutate(across(trial_offset, ~tidyr::replace_na(., mean(., na.rm=TRUE)))) %>%
  mutate(across(logRT, ~tidyr::replace_na(., mean(., na.rm=TRUE)))) %>%
  ungroup()

# Create predictor for AR-1 parameter
# 2-level model
# Create lagged outcome for AR1
jointdat = jointdat %>% 
  group_by(subject) %>%
  # neural ar-1
  mutate(Y_eeg1 = lag(trial_offset)) %>%
  mutate(Y_eeg1 = if_else(is.na(Y_eeg1), mean(Y_eeg1, na.rm=TRUE), Y_eeg1)) %>%
  # behavioral ar-1
  mutate(Y_RT1 = lag(logRT)) %>%
  mutate(Y_RT1 = if_else(is.na(Y_RT1), mean(Y_RT1, na.rm=TRUE), Y_RT1)) %>%
  ungroup()

# within-person centered time
jointdat = jointdat %>%
  group_by(subject) %>%
  mutate(ctime = trialnr/100) %>%
  mutate(mtime = mean(ctime)) %>%
  ungroup()

l2time = jointdat$mtime[which(jointdat$trialnr == "1")]

#0.4 Specify data list for stan to process
dsem_list <- list(N_subj = length(unique(as.numeric(as.factor(jointdat$subject)))), # subject number
                  Y_eeg = jointdat$trial_offset, # outcome variable EEG
                  Y_RT = jointdat$logRT, # outcome variable RT
                  N_obs = length(unique(jointdat$trialnr)),
                  subj = as.numeric(as.factor(jointdat$subject)),
                  N = nrow(jointdat),
                  time = jointdat$ctime,
                  mtime = l2time,
                  Y_eeg1 = jointdat$Y_eeg1,
                  Y_RT1 = jointdat$Y_RT1)

##############################################################################
#1. Population coding model
##############################################################################

fast_stan1=
  
  "data {
  int<lower = 1> N;
  int<lower = 1> N_subj;
  array[N] int<lower=1, upper=N_subj> subj;
  vector[N] Y_eeg;
  vector[N] Y_RT;
}

parameters {
  real<lower=0> sigma1;
  real sigma2;
  vector<lower=0>[3] tau_u;
  real alpha1;
  real alpha2;
  matrix[3, N_subj] z_u;
  cholesky_factor_corr[3] L_u;
}

transformed parameters {
matrix[N_subj,3] u;
u = transpose(diag_pre_multiply(tau_u, L_u) * z_u);
}

model {
  // priors EEG
  alpha1 ~ normal(0, 50);
  sigma1 ~ normal(0, 50);
  // priors RT
  alpha2 ~ normal(0, 50);
  sigma2 ~ normal(0, 50);
  tau_u ~ cauchy(0,2);
  to_vector(z_u) ~ std_normal();
  L_u ~ lkj_corr_cholesky(1);
  // Specify model
  Y_eeg ~ normal(alpha1 + u[subj,1], sigma1);
  Y_RT ~ normal(alpha2 + u[subj,2], exp(sigma2 + u[subj,3]));
}

generated quantities {
corr_matrix[3] rho_u = L_u * L_u';
vector[N] log_lik;
vector[N] eeg_pred;
vector[N] rt_pred;
for (i in 1:N){
log_lik[i] = normal_lpdf(Y_eeg[i] | alpha1 + u[subj[i],1], sigma1) +
normal_lpdf(Y_RT[i] | alpha2 + u[subj[i],2], exp(sigma2 + u[subj[i],3]));
}

for (i in 1:N){
eeg_pred[i] = normal_rng(alpha1 + u[subj[i],1], sigma1);
}

for (i in 1:N){
rt_pred[i] = normal_rng(alpha2 + u[subj[i],2], exp(sigma2 + u[subj[i],3]));
}


}
"

####################################
#1.1. Compile and sample from DSEM
####################################

# Save to text file
#writeLines(fast_stan1, "dsem_popcode.stan")
# Compile model
model1 <- cmdstan_model(stan_file = "dsem_popcode.stan", stanc_options = list("O1"))
# Sample from the model
#dsem_list$grainsize <- 1 #set grainsize (how many subjects per thread)
time1 = system.time(fit1 <- model1$sample(dsem_list,
                                          chains = 4,
                                          parallel_chains = 4,
                                          #threads_per_chain = 4,
                                          iter_warmup = 2000,
                                          iter_sampling = 4000,
                                          save_warmup = FALSE))

######################################
#1.2. Describe results
######################################

# Summary statistics
sum1<-fit1$summary(c("alpha1", "sigma1", "alpha2",
                     "sigma2", "tau_u", "rho_u"))
b1_draws<-fit1$draws(inc_warmup = FALSE, format = "df")
print(sum1, n = 24)

#######################################
#1.3. Model diagnostics
#######################################

# Diagnostic plots
color_scheme_set("viridisE")
trace<-mcmc_trace(b1_draws, pars = c("alpha1", "sigma1", "alpha2",
                                     "sigma2",
                                     "tau_u[1]", "tau_u[2]", "tau_u[3]"))
hist<-mcmc_hist(b1_draws, pars = c("alpha1", "sigma1", "alpha2",
                                   "sigma2",
                                   "tau_u[1]", "tau_u[2]", "tau_u[3]"))
acf<-mcmc_acf_bar(b1_draws, pars = c("alpha1", "sigma1", "alpha2",
                                     "sigma2",
                                     "tau_u[1]", "tau_u[2]", "tau_u[3]"))


# name plot files
trace_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\trace_pc_", regions, ".png")
hist_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\hist_pc_", regions, ".png")
acf_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\acf_pc_", regions, ".png")
# save plots
ggsave(plot = trace, trace_name, dpi=300, height = 10, width = 15)
ggsave(plot = hist, hist_name, dpi=300, height = 10, width = 15)
ggsave(plot = acf, acf_name, dpi=300, height = 10, width = 15)
# Free up space
rm(trace,hist,acf)
gc()
# Alert that tasks are done
beep(sound = 11, expr = NULL)

###################################################################
#1.4. Posterior predictive checks
###################################################################

eeg_pred <- as.matrix(b1_draws %>% dplyr::select(starts_with("eeg_pred")))
rt_pred <- as.matrix(b1_draws %>% dplyr::select(starts_with("rt_pred")))

# Histograms of replicated data (EEG data)
############################################
pc_hist_eeg<-ppc_hist(jointdat$trial_offset, eeg_pred[1:8, ], binwidth = 0.01)
# Compare density estimate of y to estimates from a bunch of y_rep
pc_dens_eeg<-ppc_dens_overlay(jointdat$trial_offset, eeg_pred[1:50, ])
# Scatterplot
pc_stat_eeg<-ppc_stat_2d(jointdat$trial_offset, eeg_pred, stat = c("mean", "sd"))

# name plot files
hist_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\pc_hist_eeg", regions, ".png")
dens_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\pc_dens_eeg", regions, ".png")
pcstat_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\pc_stat_eeg", regions, ".png")
# Save plots
ggsave(plot = pc_hist_eeg, hist_name, dpi=300, height = 10, width = 15)
ggsave(plot = pc_dens_eeg, dens_name, dpi=300, height = 10, width = 15)
ggsave(plot = pc_stat_eeg, pcstat_name, dpi=300, height = 10, width = 15)

rm(pc_hist_eeg, pc_dens_eeg, pc_stat_eeg)

# Histograms of replicated data (RT data)
###########################################
pc_hist_rt<-ppc_hist(jointdat$logRT, rt_pred[1:8, ], binwidth = 0.01)
# Compare density estimate of y to estimates from a bunch of y_rep
pc_dens_rt<-ppc_dens_overlay(jointdat$logRT, rt_pred[1:50, ])
# Scatterplot
pc_stat_rt<-ppc_stat_2d(jointdat$logRT, rt_pred, stat = c("mean", "sd"))

# name plot files
hist_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\pc_hist_rt", regions, ".png")
dens_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\pc_dens_rt", regions, ".png")
pcstat_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\pc_stat_rt", regions, ".png")
# Save plots
ggsave(plot = pc_hist_rt, hist_name, dpi=300, height = 10, width = 15)
ggsave(plot = pc_dens_rt, dens_name, dpi=300, height = 10, width = 15)
ggsave(plot = pc_stat_rt, pcstat_name, dpi=300, height = 10, width = 15)

rm(pc_hist_rt, pc_dens_rt, pc_stat_rt)

# Lineplot for posterior predictive check (EEG)
############################################################

#1.4.1.1 For posterior predictive data
# Estimate mean of all columns
eeg_mean = as.data.frame(colMeans(eeg_pred))
# Create a counter that increases every 100 rows
eeg_mean <- eeg_mean %>%
  mutate(subject = ceiling(row_number() / 100))
eeg_mean$eeg_est = eeg_mean$`colMeans(eeg_pred)` 
# Create a trial counter
eeg_mean<-eeg_mean %>%
  dplyr::group_by(subject) %>%
  dplyr::mutate(grp = cumsum(row_number(1))) %>%
  dplyr::group_by(subject,grp) %>%
  dplyr::mutate(trialnr = row_number()) %>%
  ungroup() %>%
  select(-grp)
# Create a trial average estimate
eeg_mean <- eeg_mean %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(trial_eeg = mean(eeg_est)) %>%
  ungroup()
# Estimate quantiles
eeg_mean <- eeg_mean %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(lower_eeg_pred = quantile(eeg_est, 0.025)) %>%
  dplyr::mutate(upper_eeg_pred = quantile(eeg_est, 0.975)) %>%
  ungroup()

#1.4.1.2 For raw data
# Create a trial average estimate
jointdat <- jointdat %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(plot_offset = mean(trial_offset)) %>%
  ungroup()

# Limit to only 1 row per trial (the estimates are the mean over all people)
plotdat_pc <- jointdat %>%
  filter(subject == 1) %>%
  select(trialnr, plot_offset)

eeg_mean = eeg_mean %>%
  filter(subject == 1) %>%
  select(trialnr, trial_eeg, lower_eeg_pred, upper_eeg_pred)

plotdat_pc <- plotdat_pc %>%
  distinct() %>%
  left_join(eeg_mean, by = "trialnr")

# Create line plot for all data
eegplot_pc = ggplot() +
  # Observed EEG (actual data)
  geom_line(data = jointdat, aes(y = trial_offset, x =  trialnr, group = as.factor(subject)), color = "black", size = 1, alpha = 0.4) +
  # Mean observed EEG (actual data)
  geom_line(data = plotdat_pc, aes(y = plot_offset, x =  trialnr), color = "#b73779", size = 1) +
  # Posterior predicted EEG (mean prediction)
  geom_line(data = plotdat_pc, aes(y = trial_eeg, x = trialnr), color = "#21918c", size = 1) +
  # Add credible intervals (95%)
  geom_ribbon(data = plotdat_pc, aes(x = trialnr, y = trial_eeg, ymin = lower_eeg_pred, ymax = upper_eeg_pred), alpha = 0.2, fill = "#21918c") +
  labs(title = "PPC for Population Coding model (EEG)",
       x = "Trial Number",
       y = "1/f Offset") +
  theme_minimal(base_size = 12)

eeg_pc <- paste0("~\\EEGvar\\Visualization\\eegplot_pc_", regions, ".png")
ggsave(plot = eegplot_pc, eeg_pc, dpi=300,
       height = 8, width = 6, bg = "white")

# Create line plot for mean only
eegplot_pc = ggplot() +
  # Observed EEG (actual data)
  #geom_line(data = jointdat, aes(y = trial_offset, x =  trialnr, group = as.factor(subject)), color = "black", size = 1, alpha = 0.4) +
  # Mean observed EEG (actual data)
  geom_line(data = plotdat_pc, aes(y = plot_offset, x =  trialnr), color = "#b73779", size = 1) +
  # Posterior predicted EEG (mean prediction)
  geom_line(data = plotdat_pc, aes(y = trial_eeg, x = trialnr), color = "#21918c", size = 1) +
  # Add credible intervals (95%)
  #geom_ribbon(data = plotdat, aes(x = trialnr, y = trial_eeg, ymin = lower_eeg_pred, ymax = upper_eeg_pred), alpha = 0.2, fill = "#fcfdbf") +
  labs(title = "PPC for Population Coding model (EEG)",
       x = "Trial Number",
       y = "1/f Offset") +
  theme_minimal(base_size = 12)+
  theme(
    plot.title = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white"),
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    panel.grid = element_line(color = "gray")
  )

eeg_pc_mean <- paste0("~\\EEGvar\\Visualization\\eegplot_pc_mean", regions, ".png")
ggsave(plot = eegplot_pc, eeg_pc_mean, dpi=300,
       height = 8, width = 6, bg = "black")

rm(eegplot_pc)

# Lineplot for posterior predictive check (RT)
################################################

#1.4.2.1 For posterior predictive data
# Estimate mean of all columns
rt_mean = as.data.frame(colMeans(rt_pred))
# Create a counter that increases every 100 rows
rt_mean <- rt_mean %>%
  mutate(subject = ceiling(row_number() / 100))
rt_mean$rt_est = rt_mean$`colMeans(rt_pred)` 
# Create a trial counter
rt_mean<-rt_mean %>%
  dplyr::group_by(subject) %>%
  dplyr::mutate(grp = cumsum(row_number(1))) %>%
  dplyr::group_by(subject,grp) %>%
  dplyr::mutate(trialnr = row_number()) %>%
  ungroup() %>%
  select(-grp)
# Create a trial average estimate
rt_mean <- rt_mean %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(trial_rt = mean(rt_est)) %>%
  ungroup()
# Estimate quantiles
rt_mean <- rt_mean %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(lower_rt_pred = quantile(rt_est, 0.025)) %>%
  dplyr::mutate(upper_rt_pred = quantile(rt_est, 0.975)) %>%
  ungroup()

#1.4.3.2 For raw data
# Create a trial average estimate
jointdat <- jointdat %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(plot_rt = mean(logRT)) %>%
  ungroup()

# Limit to only 1 row per trial (the estimates are the mean over all people)
plotdat_pcrt <- jointdat %>%
  filter(subject == 1) %>%
  select(trialnr, plot_rt)

rt_mean = rt_mean %>%
  filter(subject == 1) %>%
  select(trialnr, trial_rt, lower_rt_pred, upper_rt_pred)

plotdat_pcrt <- plotdat_pcrt %>%
  distinct() %>%
  left_join(rt_mean, by = "trialnr")

# Plot for the mean only
# Create line plot
rtplot_pc = ggplot() +
  # Observed rt (actual data)
  #geom_line(data = jointdat, aes(y = logRT, x = trialnr, group = as.factor(subject)), color = "black", size = 1, alpha = 0.4) +
  # Mean observed rt (actual data)
  geom_line(data = plotdat_pcrt, aes(y = plot_rt, x = trialnr), color = "#b73779", size = 1) +
  # Posterior predicted rt (mean prediction)
  geom_line(data = plotdat_pcrt, aes(y = trial_rt, x = trialnr), color = "#21918c", size = 1) +
  # Add credible intervals (95%)
  #geom_ribbon(data = plotdat, aes(x = trialnr, y= trial_rt, ymin = lower_rt_pred, ymax = upper_rt_pred), alpha = 0.2, fill = "#fcfdbf") +
  labs(title = "PPC for Population Coding model (RT)",
       x = "Trial Number",
       y = "Response Time") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white"),
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    panel.grid = element_line(color = "gray")
  )

rtplot = paste0("~\\EEGvar\\Visualization\\rtplot_pc_mean", regions, ".png")
ggsave(plot = rtplot_pc, rtplot, dpi=300,
       height = 8, width = 6, bg = "black")

# Plot for all participants
# Create line plot
rtplot_pc = ggplot() +
  # Observed rt (actual data)
  geom_line(data = jointdat, aes(y = logRT, x = trialnr, group = as.factor(subject)), color = "black", size = 1, alpha = 0.4) +
  # Mean observed rt (actual data)
  geom_line(data = plotdat_pcrt, aes(y = plot_rt, x = trialnr), color = "#b73779", size = 1) +
  # Posterior predicted rt (mean prediction)
  geom_line(data = plotdat_pcrt, aes(y = trial_rt, x = trialnr), color = "#21918c", size = 1) +
  # Add credible intervals (95%)
  geom_ribbon(data = plotdat_pcrt, aes(x = trialnr, y= trial_rt, ymin = lower_rt_pred, ymax = upper_rt_pred), alpha = 0.2, fill = "#21918c") +
  labs(title = "PPC for Population Coding model (RT)",
       x = "Trial Number",
       y = "Response Time") +
  theme_minimal(base_size = 12)

rtplot = paste0("~\\EEGvar\\Visualization\\rtplot_pc_alldata", regions, ".png")
ggsave(plot = rtplot_pc, rtplot, dpi=300,
       height = 8, width = 6, bg = "white")

rm(rtplot_pc)

##############################################################################
#1.5. Create correlation matrix to visualize between-subject covariance
##############################################################################

# extract between-subject correlation matrix
pc_covmat = as.matrix(b1_draws %>% dplyr::select(starts_with("rho_u")))
# take the mean of correlations
pc_meancov = apply(pc_covmat, MARGIN = c(2), mean)
pc_meancov = matrix(pc_meancov, nrow = 3, ncol = 3, byrow = TRUE)
rownames(pc_meancov) = colnames(pc_meancov) = c("alpha1", "alpha2", "psi2")
# plot correlation matrix
# Set up the PNG file
corrname = paste0("~\\EEGvar\\Visualization\\pc_corrplot", regions,".png")
png(corrname, res = 300, width = 1000, height = 1000)  # Specify file name and size

corrplot(pc_meancov, method="color", 
         addCoef.col = "black",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         type = "lower",
         diag = TRUE,
         tl.col = "black")

dev.off()
gc()

##############################################################################
#2. Variable precision model
##############################################################################

fast_stan2=
  
  "data {
  int<lower = 1> N;
  int<lower = 1> N_subj;
  array[N] int<lower=1, upper=N_subj> subj;
  vector[N] Y_eeg;
  vector[N] Y_RT;
}

parameters {
  real sigma1;
  real sigma2;
  vector<lower=0>[5] tau_u;
  real alpha1;
  real alpha2;
  real couple3;
  matrix[5, N_subj] z_u;
  cholesky_factor_corr[5] L_u;
}

transformed parameters {
matrix[N_subj,5] u;
u = transpose(diag_pre_multiply(tau_u, L_u) * z_u);
}

model {
  // priors EEG
  alpha1 ~ normal(0, 50);
  sigma1 ~ normal(0, 50);
  // priors RT
  alpha2 ~ normal(0, 50);
  sigma2 ~ normal(0, 50);
  couple3 ~ normal(0,50);
  tau_u ~ cauchy(0,2);
  to_vector(z_u) ~ std_normal();
  L_u ~ lkj_corr_cholesky(1);
  // Specify model
  Y_eeg ~ normal(alpha1 + u[subj,1], exp(sigma1 + u[subj,2]));
  Y_RT ~ normal(alpha2 + u[subj,3] + (Y_eeg - (alpha1 + u[subj,1])).*(couple3+u[subj,4]), 
              exp(sigma2 + u[subj,5]));
}

generated quantities {
corr_matrix[5] rho_u = L_u * L_u';
vector[N] log_lik;
vector[N] eeg_pred;
vector[N] rt_pred;
for (i in 1:N){
log_lik[i] = normal_lpdf(Y_eeg[i] | alpha1 + u[subj[i],1], exp(sigma1 + u[subj[i],2])) +
normal_lpdf(Y_RT[i] | alpha2 + u[subj[i],3] + (Y_eeg[i] - (alpha1 + u[subj[i],1])).*(couple3+u[subj[i],4]), 
              exp(sigma2 + u[subj[i],5]));
}

for (i in 1:N){
eeg_pred[i] = normal_rng(alpha1 + u[subj[i],1], exp(sigma1 + u[subj[i],2]));
}

for (i in 1:N){
rt_pred[i] = normal_rng(alpha2 + u[subj[i],3] + (Y_eeg[i] - (alpha1 + u[subj[i],1])).*(couple3+u[subj[i],4]), 
              exp(sigma2 + u[subj[i],5]));
}

}
"

####################################
#2.1. Compile and sample from DSEM
####################################

# Save to text file
#writeLines(fast_stan2, "dsem_vpmodel.stan")
# Compile model
model2 <- cmdstan_model(stan_file = "dsem_vpmodel.stan", stanc_options = list("O1"))
# Sample from the model
#dsem_list$grainsize <- 1 #set grainsize (how many subjects per thread)
time2 = system.time(fit2 <- model2$sample(dsem_list,
                                          chains = 4,
                                          parallel_chains = 4,
                                          #threads_per_chain = 4,
                                          iter_warmup = 2000,
                                          iter_sampling = 4000))

######################################
#2.2. Describe results
######################################

# Summary statistics
sum2<-fit2$summary(c("alpha1", "sigma1", "alpha2",
                     "sigma2", "couple3",
                     "tau_u", "rho_u"))
b2_draws<-fit2$draws(format = "df")
print(sum2, n = 24)

#######################################
#2.3. Model diagnostics
#######################################

# Diagnostic plots
color_scheme_set("viridisE")
trace2<-mcmc_trace(b2_draws, pars = c("alpha1", "sigma1", "alpha2",
                                     "sigma2", "couple3",
                                     "tau_u[1]", "tau_u[2]", "tau_u[3]",
                                     "tau_u[4]", "tau_u[5]"))
hist2<-mcmc_hist(b2_draws, pars = c("alpha1", "sigma1", "alpha2",
                                   "sigma2", "couple3",
                                   "tau_u[1]", "tau_u[2]", "tau_u[3]",
                                   "tau_u[4]", "tau_u[5]"))
acf2<-mcmc_acf_bar(b2_draws, pars = c("alpha1", "sigma1", "alpha2",
                                     "sigma2", "couple3",
                                     "tau_u[1]", "tau_u[2]", "tau_u[3]",
                                     "tau_u[4]", "tau_u[5]"))

# name plot files
trace_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\trace_vp_", regions, ".png")
hist_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\hist_vp_", regions, ".png")
acf_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\acf_vp_", regions, ".png")
# save plots
ggsave(plot = trace2, trace_name, dpi=300, height = 10, width = 15)
ggsave(plot = hist2, hist_name, dpi=300, height = 10, width = 15)
ggsave(plot = acf2, acf_name, dpi=300, height = 10, width = 15)

# Free up space
rm(trace2,hist2,acf2)
gc()
beep(sound = 11, expr = NULL)

###################################################################
#2.4. Posterior predictive checks
###################################################################

eeg_pred <- as.matrix(b2_draws %>% dplyr::select(starts_with("eeg_pred")))
rt_pred <- as.matrix(b2_draws %>% dplyr::select(starts_with("rt_pred")))

# Histograms of replicated data (EEG data)
############################################
vp_hist_eeg<-ppc_hist(jointdat$trial_offset, eeg_pred[1:8, ], binwidth = 0.01)
# Compare density estimate of y to estimates from a bunch of y_rep
vp_dens_eeg<-ppc_dens_overlay(jointdat$trial_offset, eeg_pred[1:50, ])
# Scatter plot
vp_stat_eeg<-ppc_stat_2d(jointdat$trial_offset, eeg_pred, stat = c("mean", "sd"))

# name plot files
hist_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\vp_hist_eeg", regions, ".png")
dens_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\vp_dens_eeg", regions, ".png")
pcstat_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\vp_stat_eeg", regions, ".png")
# Save plots
ggsave(plot = vp_hist_eeg, hist_name, dpi=300, height = 10, width = 15)
ggsave(plot = vp_dens_eeg, dens_name, dpi=300, height = 10, width = 15)
ggsave(plot = vp_stat_eeg, pcstat_name, dpi=300, height = 10, width = 15)

rm(vp_hist_eeg, vp_dens_eeg, vp_stat_eeg)

# Histograms of replicated data (RT data)
###########################################
vp_hist_rt<-ppc_hist(jointdat$logRT, rt_pred[1:8, ], binwidth = 0.01)
# Compare density estimate of y to estimates from a bunch of y_rep
vp_dens_rt<-ppc_dens_overlay(jointdat$logRT, rt_pred[1:50, ])
# Scatterplot
vp_stat_rt<-ppc_stat_2d(jointdat$logRT, rt_pred, stat = c("mean", "sd"))

# name plot files
hist_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\vp_hist_rt", regions, ".png")
dens_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\vp_dens_rt", regions, ".png")
pcstat_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\vp_stat_rt", regions, ".png")
# Save plots
ggsave(plot = vp_hist_rt, hist_name, dpi=300, height = 10, width = 15)
ggsave(plot = vp_dens_rt, dens_name, dpi=300, height = 10, width = 15)
ggsave(plot = vp_stat_rt, pcstat_name, dpi=300, height = 10, width = 15)

rm(vp_hist_rt, vp_dens_rt, vp_stat_rt)

# Lineplot for posterior predictive check (EEG)
############################################################

#1.4.1.1 For posterior predictive data
# Estimate mean of all columns
eeg_mean = as.data.frame(colMeans(eeg_pred))
# Create a counter that increases every 100 rows
eeg_mean <- eeg_mean %>%
  mutate(subject = ceiling(row_number() / 100))
eeg_mean$eeg_est = eeg_mean$`colMeans(eeg_pred)` 
# Create a trial counter
eeg_mean<-eeg_mean %>%
  dplyr::group_by(subject) %>%
  dplyr::mutate(grp = cumsum(row_number(1))) %>%
  dplyr::group_by(subject,grp) %>%
  dplyr::mutate(trialnr = row_number()) %>%
  ungroup() %>%
  select(-grp)
# Create a trial average estimate
eeg_mean <- eeg_mean %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(trial_eeg = mean(eeg_est)) %>%
  ungroup()
# Estimate quantiles
eeg_mean <- eeg_mean %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(lower_eeg_pred = quantile(eeg_est, 0.025)) %>%
  dplyr::mutate(upper_eeg_pred = quantile(eeg_est, 0.975)) %>%
  ungroup()

#1.4.1.2 For raw data
# Create a trial average estimate
jointdat <- jointdat %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(plot_offset = mean(trial_offset)) %>%
  ungroup()

# Limit to only 1 row per trial (the estimates are the mean over all people)
plotdat_vp <- jointdat %>%
  filter(subject == 1) %>%
  select(trialnr, plot_offset)

eeg_mean = eeg_mean %>%
  filter(subject == 1) %>%
  select(trialnr, trial_eeg, lower_eeg_pred, upper_eeg_pred)

plotdat_vp <- plotdat_vp %>%
  distinct() %>%
  left_join(eeg_mean, by = "trialnr")

# Create line plot for all data
eegplot_vp = ggplot() +
  # Observed EEG (actual data)
  geom_line(data = jointdat, aes(y = trial_offset, x =  trialnr, group = as.factor(subject)), color = "black", size = 1, alpha = 0.4) +
  # Mean observed EEG (actual data)
  geom_line(data = plotdat_vp, aes(y = plot_offset, x =  trialnr), color = "#b73779", size = 1) +
  # Posterior predicted EEG (mean prediction)
  geom_line(data = plotdat_vp, aes(y = trial_eeg, x = trialnr), color = "#fc8961", size = 1) +
  # Add credible intervals (95%)
  geom_ribbon(data = plotdat_vp, aes(x = trialnr, y = trial_eeg, ymin = lower_eeg_pred, ymax = upper_eeg_pred), alpha = 0.2, fill = "#fc8961") +
  labs(title = "PPC for Variable Precision model (EEG)",
       x = "Trial Number",
       y = "1/f Offset") +
  theme_minimal(base_size = 12)

eeg_vp <- paste0("~\\EEGvar\\Visualization\\eegplot_vp_alldat", regions, ".png")
ggsave(plot = eegplot_vp, eeg_vp, dpi=300,
       height = 8, width = 6, bg = "white")

# Create line plot for mean only
eegplot_vp = ggplot() +
  # Observed EEG (actual data)
  #geom_line(data = jointdat, aes(y = trial_offset, x =  trialnr, group = as.factor(subject)), color = "black", size = 1, alpha = 0.4) +
  # Mean observed EEG (actual data)
  geom_line(data = plotdat_vp, aes(y = plot_offset, x =  trialnr), color = "#b73779", size = 1) +
  # Posterior predicted EEG (mean prediction)
  geom_line(data = plotdat_vp, aes(y = trial_eeg, x = trialnr), color = "#fc8961", size = 1) +
  # Add credible intervals (95%)
  #geom_ribbon(data = plotdat, aes(x = trialnr, y = trial_eeg, ymin = lower_eeg_pred, ymax = upper_eeg_pred), alpha = 0.2, fill = "#fcfdbf") +
  labs(title = "PPC for Variable Precision model (EEG)",
       x = "Trial Number",
       y = "1/f Offset") +
  theme_minimal(base_size = 12)+
  theme(
    plot.title = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white"),
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    panel.grid = element_line(color = "gray")
  )

eeg_vp <- paste0("~\\EEGvar\\Visualization\\eegplot_vp_mean", regions, ".png")
ggsave(plot = eegplot_vp, eeg_vp, dpi=300,
       height = 8, width = 6, bg = "black")

rm(eegplot_vp)

# Lineplot for posterior predictive check (RT)
################################################

#1.4.2.1 For posterior predictive data
# Estimate mean of all columns
rt_mean = as.data.frame(colMeans(rt_pred))
# Create a counter that increases every 100 rows
rt_mean <- rt_mean %>%
  mutate(subject = ceiling(row_number() / 100))
rt_mean$rt_est = rt_mean$`colMeans(rt_pred)` 
# Create a trial counter
rt_mean<-rt_mean %>%
  dplyr::group_by(subject) %>%
  dplyr::mutate(grp = cumsum(row_number(1))) %>%
  dplyr::group_by(subject,grp) %>%
  dplyr::mutate(trialnr = row_number()) %>%
  ungroup() %>%
  select(-grp)
# Create a trial average estimate
rt_mean <- rt_mean %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(trial_rt = mean(rt_est)) %>%
  ungroup()
# Estimate quantiles
rt_mean <- rt_mean %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(lower_rt_pred = quantile(rt_est, 0.025)) %>%
  dplyr::mutate(upper_rt_pred = quantile(rt_est, 0.975)) %>%
  ungroup()

#1.4.3.2 For raw data
# Create a trial average estimate
jointdat <- jointdat %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(plot_rt = mean(logRT)) %>%
  ungroup()

# Limit to only 1 row per trial (the estimates are the mean over all people)
plotdat_vprt <- jointdat %>%
  filter(subject == 1) %>%
  select(trialnr, plot_rt)

rt_mean = rt_mean %>%
  filter(subject == 1) %>%
  select(trialnr, trial_rt, lower_rt_pred, upper_rt_pred)

plotdat_vprt <- plotdat_vprt %>%
  distinct() %>%
  left_join(rt_mean, by = "trialnr")

# Plot for the mean only
# Create line plot
rtplot_vp = ggplot() +
  # Observed rt (actual data)
  #geom_line(data = jointdat, aes(y = logRT, x = trialnr, group = as.factor(subject)), color = "black", size = 1, alpha = 0.4) +
  # Mean observed rt (actual data)
  geom_line(data = plotdat_vprt, aes(y = plot_rt, x = trialnr), color = "#b73779", size = 1) +
  # Posterior predicted rt (mean prediction)
  geom_line(data = plotdat_vprt, aes(y = trial_rt, x = trialnr), color = "#fc8961", size = 1) +
  # Add credible intervals (95%)
  #geom_ribbon(data = plotdat, aes(x = trialnr, y= trial_rt, ymin = lower_rt_pred, ymax = upper_rt_pred), alpha = 0.2, fill = "#fcfdbf") +
  labs(title = "PPC for Variable Precision model (RT)",
       x = "Trial Number",
       y = "Response Time") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white"),
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    panel.grid = element_line(color = "gray")
  )

rt_vp = paste0("~\\EEGvar\\Visualization\\rtplot_vp_mean", regions, ".png")
ggsave(plot = rtplot_vp, rt_vp, dpi=300,
       height = 8, width = 6, bg = "black")

# Plot for all participants
# Create line plot
rtplot_vp = ggplot() +
  # Observed rt (actual data)
  geom_line(data = jointdat, aes(y = logRT, x = trialnr, group = as.factor(subject)), color = "black", size = 1, alpha = 0.4) +
  # Mean observed rt (actual data)
  geom_line(data = plotdat_vprt, aes(y = plot_rt, x = trialnr), color = "#b73779", size = 1) +
  # Posterior predicted rt (mean prediction)
  geom_line(data = plotdat_vprt, aes(y = trial_rt, x = trialnr), color = "#fc8961", size = 1) +
  # Add credible intervals (95%)
  geom_ribbon(data = plotdat_vprt, aes(x = trialnr, y= trial_rt, ymin = lower_rt_pred, ymax = upper_rt_pred), alpha = 0.2, fill = "#fc8961") +
  labs(title = "PPC for Variable Precision model (RT)",
       x = "Trial Number",
       y = "Response Time") +
  theme_minimal(base_size = 12)

rt_vp = paste0("~\\EEGvar\\Visualization\\rtplot_vp_alldata", regions, ".png")
ggsave(plot = rtplot_vp, rt_vp, dpi=300,
       height = 8, width = 6, bg = "white")

rm(rtplot_vp)

##############################################################################
#2.5. Create correlation matrix to visualize between-subject covariance
##############################################################################

# extract between-subject correlation matrix
vp_covmat = as.matrix(b2_draws %>% dplyr::select(starts_with("rho_u")))
# take the mean of correlations
vp_meancov = apply(vp_covmat, MARGIN = c(2), mean)
vp_meancov = matrix(vp_meancov, nrow = 5, ncol = 5, byrow = TRUE)
rownames(vp_meancov) = colnames(vp_meancov) = c("alpha1", "psi1", "alpha2",
                                                "gamma1", "psi2")
# plot correlation matrix
# Set up the PNG file
corrfig = paste0("~\\EEGvar\\Visualization\\vp_corrplot", regions, ".png")
png(corrfig, res = 300, width = 1000, height = 1000)  # Specify file name and size

corrplot(vp_meancov, method="color", 
         addCoef.col = "black",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         type = "lower",
         diag = TRUE,
         tl.col = "black")

dev.off()
gc()

##############################################################################
#3. Sequence sensitive model
##############################################################################

fast_stan3=
  
  "data {
  int<lower = 1> N;
  int<lower = 1> N_subj;
  array[N] int<lower=1, upper=N_subj> subj;
  vector[N] Y_eeg;
  vector[N] Y_RT;
  vector[N] time;
  vector[N] Y_eeg1;
  vector[N] Y_RT1;
  vector[N_subj] mtime;
}

parameters {
  real sigma1;
  real sigma2;
  vector<lower=0>[10] tau_u;
  real alpha1;
  real beta1;
  real phi1;
  real alpha2;
  real beta2;
  real phi2;
  real couple1;
  real couple3;
  matrix[10, N_subj] z_u;
  cholesky_factor_corr[10] L_u;
}

transformed parameters {
matrix[N_subj,10] u;
u = transpose(diag_pre_multiply(tau_u, L_u) * z_u);
}

model {
  // priors EEG
  alpha1 ~ normal(0, 50);
  beta1 ~ normal(0,50);
  phi1 ~ normal(0,50);
  sigma1 ~ normal(0, 50);
  couple1 ~ normal(0,50);
  // priors RT
  alpha2 ~ normal(0, 50);
  phi2 ~ normal(0,50);
  beta2 ~ normal(0,50);
  sigma2 ~ normal(0, 50);
  couple3 ~ normal(0,50);
  tau_u ~ cauchy(0,2);
  to_vector(z_u) ~ std_normal();
  L_u ~ lkj_corr_cholesky(1);
  // Specify model
  Y_eeg ~ normal(alpha1 + u[subj,1] + (time - mtime[subj]).*(beta1+u[subj,2]) + (Y_eeg1 - (alpha1 + u[subj,1])).*(phi1+u[subj,3]) + (Y_RT1 - (alpha2 + u[subj,6])).*(couple1+u[subj,4]), 
              exp(sigma1 + u[subj,5]));
  Y_RT ~ normal(alpha2 + u[subj,6] + (time - mtime[subj]).*(beta2+u[subj,7]) + (Y_RT1 - (alpha2 + u[subj,6])).*(phi2+u[subj,8]) + (Y_eeg - (alpha1 + u[subj,1])).*(couple3+u[subj,9]), 
              exp(sigma2 + u[subj,10]));
}

generated quantities {
corr_matrix[10] rho_u = L_u * L_u';
vector[N] log_lik;
vector[N] eeg_pred;
vector[N] rt_pred;
for (i in 1:N){
  log_lik[i] = normal_lpdf(Y_eeg[i] | alpha1 + u[subj[i],1] + (time[i] - mtime[subj[i]]).*(beta1+u[subj[i],2]) + (Y_eeg1[i] - (alpha1 + u[subj[i],1])).*(phi1+u[subj[i],3]) + (Y_RT1[i] - (alpha2 + u[subj[i],6])).*(couple1+u[subj[i],4]), 
              exp(sigma1 + u[subj[i],5])) +
  normal_lpdf(Y_RT[i] | alpha2 + u[subj[i],6] + (time[i] - mtime[subj[i]]).*(beta2+u[subj[i],7]) + (Y_RT1[i] - (alpha2 + u[subj[i],6])).*(phi2+u[subj[i],8]) + (Y_eeg[i] - (alpha1 + u[subj[i],1])).*(couple3+u[subj[i],9]),
              exp(sigma2 + u[subj[i],10]));
}

for (i in 1:N){
eeg_pred[i] = normal_rng(alpha1 + u[subj[i],1] + (time[i] - mtime[subj[i]]).*(beta1+u[subj[i],2]) + (Y_eeg1[i] - (alpha1 + u[subj[i],1])).*(phi1+u[subj[i],3]) + (Y_RT1[i] - (alpha2 + u[subj[i],6])).*(couple1+u[subj[i],4]), 
              exp(sigma1 + u[subj[i],5]));
}

for (i in 1:N){
rt_pred[i] = normal_rng(alpha2 + u[subj[i],6] + (time[i] - mtime[subj[i]]).*(beta2+u[subj[i],7]) + (Y_RT1[i] - (alpha2 + u[subj[i],6])).*(phi2+u[subj[i],8]) + (Y_eeg[i] - (alpha1 + u[subj[i],1])).*(couple3+u[subj[i],9]),
              exp(sigma2 + u[subj[i],10]));
}

}
"

####################################
#3.1. Compile and sample from DSEM
####################################

# Save to text file
# writeLines(fast_stan3, "dsem_maximal.stan")
# Compile model
model3 <- cmdstan_model(stan_file = "dsem_maximal.stan", stanc_options = list("O1"))
# Sample from the model
#dsem_list$grainsize <- 1 #set grainsize (how many subjects per thread)
time3 = system.time(fit3 <- model3$sample(dsem_list,
                                          chains = 4,
                                          parallel_chains = 4,
                                          #threads_per_chain = 4,
                                          iter_warmup = 2000,
                                          iter_sampling = 4000))

######################################
#3.2. Describe results
######################################

# Summary statistics
sum3<-fit3$summary(c("alpha1", "beta1", "phi1", "sigma1",
                     "alpha2", "beta2", "phi2", "sigma2",
                     "couple1", "couple3",
                     "tau_u", "rho_u"))
b3_draws<-fit3$draws(format = "df")
print(sum3, n = 24)

#######################################
#3.3. Model diagnostics
#######################################

# Diagnostic plots
color_scheme_set("viridisE")
trace3<-mcmc_trace(b3_draws, pars = c("alpha1", "beta1", "phi1", "sigma1", "couple1",
                                      "alpha2", "beta2", "phi2", "sigma2", "couple3",
                                     "tau_u[1]", "tau_u[2]", "tau_u[3]", "tau_u[4]",
                                     "tau_u[5]", "tau_u[6]", "tau_u[7]", "tau_u[8]",
                                     "tau_u[9]", "tau_u[10]"))
hist3<-mcmc_hist(b3_draws, pars = c("alpha1", "beta1", "phi1", "sigma1", "couple1",
                                    "alpha2", "beta2", "phi2", "sigma2", "couple3",
                                    "tau_u[1]", "tau_u[2]", "tau_u[3]", "tau_u[4]",
                                    "tau_u[5]", "tau_u[6]", "tau_u[7]", "tau_u[8]",
                                    "tau_u[9]", "tau_u[10]"))
acf3<-mcmc_acf_bar(b3_draws, pars = c("alpha1", "beta1", "phi1", "sigma1", "couple1",
                                      "alpha2", "beta2", "phi2", "sigma2", "couple3",
                                      "tau_u[1]", "tau_u[2]", "tau_u[3]", "tau_u[4]",
                                      "tau_u[5]", "tau_u[6]", "tau_u[7]", "tau_u[8]",
                                     "tau_u[9]", "tau_u[10]"))

# name plot files
trace_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\trace_ss_", regions, ".png")
hist_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\hist_ss_", regions, ".png")
acf_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\acf_ss_", regions, ".png")
# save plots
ggsave(plot = trace3, trace_name, dpi=300, height = 10, width = 15)
ggsave(plot = hist3, hist_name, dpi=300, height = 10, width = 15)
ggsave(plot = acf3, acf_name, dpi=300, height = 10, width = 15)

# Free up space
rm(hist3, trace3, acf3)
gc()
beep(sound = 11, expr = NULL)

###################################################################
#3.4. Posterior predictive checks (code from: DYNASTI)
###################################################################

eeg_pred <- as.matrix(b3_draws %>% dplyr::select(starts_with("eeg_pred")))
rt_pred <- as.matrix(b3_draws %>% dplyr::select(starts_with("rt_pred")))

# Histograms of replicated data (EEG data)
############################################
ss_hist_eeg<-ppc_hist(jointdat$trial_offset, eeg_pred[1:8, ], binwidth = 0.01)
# Compare density estimate of y to estimates from a bunch of y_rep
ss_dens_eeg<-ppc_dens_overlay(jointdat$trial_offset, eeg_pred[1:50, ])
# Scatterplot
ss_stat_eeg<-ppc_stat_2d(jointdat$trial_offset, eeg_pred, stat = c("mean", "sd"))

# name plot files
hist_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\ss_hist_eeg_", regions, ".png")
dens_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\ss_dens_eeg_", regions, ".png")
stat_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\ss_stat_eeg_", regions, ".png")
# Save plots
ggsave(plot = ss_hist_eeg, hist_name, dpi=300, height = 6, width = 6)
ggsave(plot = ss_dens_eeg, dens_name, dpi=300, height = 6, width = 6)
ggsave(plot = ss_stat_eeg, stat_name, dpi=300, height = 6, width = 6)

rm(ss_hist_eeg, ss_dens_eeg, ss_stat_eeg)

# Histograms of replicated data (RT data)
###########################################
ss_hist_rt<-ppc_hist(jointdat$logRT, rt_pred[1:8, ], binwidth = 0.01)
# Compare density estimate of y to estimates from a bunch of y_rep
ss_dens_rt<-ppc_dens_overlay(jointdat$logRT, rt_pred[1:50, ])
# Scatterplot
ss_stat_rt<-ppc_stat_2d(jointdat$logRT, rt_pred, stat = c("mean", "sd"))

# name plot files
hist_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\ss_hist_rt_", regions, ".png")
dens_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\ss_dens_rt_", regions, ".png")
stat_name <- paste0("~\\EEGvar\\Visualization\\Diagnostics\\ss_stat_rt_", regions, ".png")
# Save plots
ggsave(plot = ss_hist_rt, hist_name, dpi=300, height = 6, width = 6)
ggsave(plot = ss_dens_rt, dens_name, dpi=300, height = 6, width = 6)
ggsave(plot = ss_stat_rt, stat_name, dpi=300, height = 6, width = 6)

rm(ss_hist_rt, ss_dens_rt, ss_stat_rt)

# Lineplot for posterior predictive check (EEG)
############################################################

#1.4.1.1 For posterior predictive data
# Estimate mean of all columns
eeg_mean = as.data.frame(colMeans(eeg_pred))
# Create a counter that increases every 100 rows
eeg_mean <- eeg_mean %>%
  mutate(subject = ceiling(row_number() / 100))
eeg_mean$eeg_est = eeg_mean$`colMeans(eeg_pred)` 
# Create a trial counter
eeg_mean<-eeg_mean %>%
  dplyr::group_by(subject) %>%
  dplyr::mutate(grp = cumsum(row_number(1))) %>%
  dplyr::group_by(subject,grp) %>%
  dplyr::mutate(trialnr = row_number()) %>%
  ungroup() %>%
  select(-grp)
# Create a trial average estimate
eeg_mean <- eeg_mean %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(trial_eeg = mean(eeg_est)) %>%
  ungroup()
# Estimate quantiles
eeg_mean <- eeg_mean %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(lower_eeg_pred = quantile(eeg_est, 0.025)) %>%
  dplyr::mutate(upper_eeg_pred = quantile(eeg_est, 0.975)) %>%
  ungroup()

#1.4.1.2 For raw data
# Create a trial average estimate
jointdat <- jointdat %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(plot_offset = mean(trial_offset)) %>%
  ungroup()

# Limit to only 1 row per trial (the estimates are the mean over all people)
plotdat_ss <- jointdat %>%
  filter(subject == 1) %>%
  select(trialnr, plot_offset)

eeg_mean = eeg_mean %>%
  filter(subject == 1) %>%
  select(trialnr, trial_eeg, lower_eeg_pred, upper_eeg_pred)

plotdat_ss <- plotdat_ss %>%
  distinct() %>%
  left_join(eeg_mean, by = "trialnr")

# Create line plot for all data
eegplot_ss = ggplot() +
  # Observed EEG (actual data)
  geom_line(data = jointdat, aes(y = trial_offset, x =  trialnr, group = as.factor(subject)), color = "black", size = 1, alpha = 0.4) +
  # Mean observed EEG (actual data)
  geom_line(data = plotdat_ss, aes(y = plot_offset, x =  trialnr), color = "#b73779", size = 1) +
  # Posterior predicted EEG (mean prediction)
  geom_line(data = plotdat_ss, aes(y = trial_eeg, x = trialnr), color = "#fcfdbf", size = 1) +
  # Add credible intervals (95%)
  geom_ribbon(data = plotdat_ss, aes(x = trialnr, y = trial_eeg, ymin = lower_eeg_pred, ymax = upper_eeg_pred), alpha = 0.2, fill = "#fcfdbf") +
  labs(title = "PPC for Sequence Sensitive model (EEG)",
       x = "Trial Number",
       y = "1/f Offset") +
  theme_minimal(base_size = 12)

eegplot = paste0("~\\EEGvar\\Visualization\\eegplot_ss_alldat_", regions, ".png")
ggsave(plot = eegplot_ss, eegplot, dpi=300,
       height = 8, width = 6, bg = "white")

# Create line plot for mean only
eegplot_ss = ggplot() +
  # Observed EEG (actual data)
  #geom_line(data = jointdat, aes(y = trial_offset, x =  trialnr, group = as.factor(subject)), color = "black", size = 1, alpha = 0.4) +
  # Mean observed EEG (actual data)
  geom_line(data = plotdat_ss, aes(y = plot_offset, x =  trialnr), color = "#b73779", size = 1) +
  # Posterior predicted EEG (mean prediction)
  geom_line(data = plotdat_ss, aes(y = trial_eeg, x = trialnr), color = "#fcfdbf", size = 1) +
  # Add credible intervals (95%)
  #geom_ribbon(data = plotdat, aes(x = trialnr, y = trial_eeg, ymin = lower_eeg_pred, ymax = upper_eeg_pred), alpha = 0.2, fill = "#fcfdbf") +
  labs(title = "PPC for Sequence Sensitive model (EEG)",
       x = "Trial Number",
       y = "1/f Offset") +
  theme_minimal(base_size = 12)+
  theme(
    plot.title = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white"),
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    panel.grid = element_line(color = "gray")
  )

eegplot = paste0("~\\EEGvar\\Visualization\\eegplot_ss_mean", regions, ".png")
ggsave(plot = eegplot_ss, eegplot, dpi=300,
       height = 8, width = 6, bg = "black")

rm(eegplot_ss)

# Lineplot for posterior predictive check (RT)
################################################

#1.4.2.1 For posterior predictive data
# Estimate mean of all columns
rt_mean = as.data.frame(colMeans(rt_pred))
# Create a counter that increases every 100 rows
rt_mean <- rt_mean %>%
  mutate(subject = ceiling(row_number() / 100))
rt_mean$rt_est = rt_mean$`colMeans(rt_pred)` 
# Create a trial counter
rt_mean<-rt_mean %>%
  dplyr::group_by(subject) %>%
  dplyr::mutate(grp = cumsum(row_number(1))) %>%
  dplyr::group_by(subject,grp) %>%
  dplyr::mutate(trialnr = row_number()) %>%
  ungroup() %>%
  select(-grp)
# Create a trial average estimate
rt_mean <- rt_mean %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(trial_rt = mean(rt_est)) %>%
  ungroup()
# Estimate quantiles
rt_mean <- rt_mean %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(lower_rt_pred = quantile(rt_est, 0.025)) %>%
  dplyr::mutate(upper_rt_pred = quantile(rt_est, 0.975)) %>%
  ungroup()

#1.4.3.2 For raw data
# Create a trial average estimate
jointdat <- jointdat %>%
  dplyr::group_by(trialnr) %>%
  dplyr::mutate(plot_rt = mean(logRT)) %>%
  ungroup()

# Limit to only 1 row per trial (the estimates are the mean over all people)
plotdat_ssrt <- jointdat %>%
  filter(subject == 1) %>%
  select(trialnr, plot_rt)

rt_mean = rt_mean %>%
  filter(subject == 1) %>%
  select(trialnr, trial_rt, lower_rt_pred, upper_rt_pred)

plotdat_ssrt <- plotdat_ssrt %>%
  distinct() %>%
  left_join(rt_mean, by = "trialnr")

# Plot for the mean only
# Create line plot
rtplot_ss = ggplot() +
  # Observed rt (actual data)
  #geom_line(data = jointdat, aes(y = logRT, x = trialnr, group = as.factor(subject)), color = "black", size = 1, alpha = 0.4) +
  # Mean observed rt (actual data)
  geom_line(data = plotdat_ssrt, aes(y = plot_rt, x = trialnr), color = "#b73779", size = 1) +
  # Posterior predicted rt (mean prediction)
  geom_line(data = plotdat_ssrt, aes(y = trial_rt, x = trialnr), color = "#fcfdbf", size = 1) +
  # Add credible intervals (95%)
  #geom_ribbon(data = plotdat, aes(x = trialnr, y= trial_rt, ymin = lower_rt_pred, ymax = upper_rt_pred), alpha = 0.2, fill = "#fcfdbf") +
  labs(title = "PPC for Sequence Sensitive model (RT)",
       x = "Trial Number",
       y = "Response Time") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white"),
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    panel.grid = element_line(color = "gray")
  )

rtplot = paste0("~\\EEGvar\\Visualization\\rtplot_ss_mean_", regions, ".png")
ggsave(plot = rtplot_ss, rtplot, dpi=300,
       height = 8, width = 6, bg = "black")

# Plot for all participants
# Create line plot
rtplot_ss = ggplot() +
  # Observed rt (actual data)
  geom_line(data = jointdat, aes(y = logRT, x = trialnr, group = as.factor(subject)), color = "black", size = 1, alpha = 0.4) +
  # Mean observed rt (actual data)
  geom_line(data = plotdat_ssrt, aes(y = plot_rt, x = trialnr), color = "#b73779", size = 1) +
  # Posterior predicted rt (mean prediction)
  geom_line(data = plotdat_ssrt, aes(y = trial_rt, x = trialnr), color = "#fcfdbf", size = 1) +
  # Add credible intervals (95%)
  geom_ribbon(data = plotdat_ssrt, aes(x = trialnr, y= trial_rt, ymin = lower_rt_pred, ymax = upper_rt_pred), alpha = 0.2, fill = "#fcfdbf") +
  labs(title = "PPC for Sequence Sensitive model (RT)",
       x = "Trial Number",
       y = "Response Time") +
  theme_minimal(base_size = 12)

rtplot = paste0("~\\EEGvar\\Visualization\\rtplot_ss_alldata_", regions, ".png")
ggsave(plot = rtplot_ss, rtplot, dpi=300,
       height = 8, width = 6, bg = "white")

rm(rtplot_ss)

##############################################################################
#3.5. Create correlation matrix to visualize between-subject covariance
##############################################################################

# extract between-subject correlation matrix
ss_covmat = as.matrix(b3_draws %>% dplyr::select(starts_with("rho_u")))
# take the mean of correlations
ss_meancov = apply(ss_covmat, MARGIN = c(2), mean)
# create nrow x ncolumn matrix
ss_meancov = matrix(ss_meancov, nrow = 10, ncol = 10, byrow = TRUE)
# write row and column names
rownames(ss_meancov) = colnames(ss_meancov) = c("alpha1", "beta1", "phi1",
                                                "gamma1", "psi1", "alpha2",
                                                "beta2", "phi2", "gamma2",
                                                "psi2")
testRes = cor.mtest(ss_meancov, conf.level = 0.95)
# plot correlation matrix
# Set up the PNG file
#cor_get_pval()
corrfig = paste0("~\\EEGvar\\Visualization\\ss_corrplot_", regions, ".png")
png(corrfig, res = 300, width = 1800, height = 1800)  # Specify file name and size
# create correlation plot
corrplot(ss_meancov, method="color", 
         addCoef.col = "black",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         type = "lower",
         diag = TRUE,
         tl.col = "black")
         #p.mat = testRes$p)

dev.off()
gc()

##################################################
#4. Model comparison (all three models)
##################################################
# Population Coding model loo
b1_ll <- as.matrix(b1_draws %>% dplyr::select(starts_with("log_lik")))
loo1<-loo::elpd(b1_ll)

# Variable Precision model loo
b2_ll <- as.matrix(b2_draws %>% dplyr::select(starts_with("log_lik")))
loo2<-loo::elpd(b2_ll)

# Sequence Sensitive model loo
b3_ll <- as.matrix(b3_draws %>% dplyr::select(starts_with("log_lik")))
loo3<-loo::elpd(b3_ll)

# Model comparison
loo::loo_compare(loo1,loo2,loo3) # compare elpds

rm(b1_ll, b2_ll, b3_ll)
gc()

##################################################
#5. Overlay all three PPC (mean lineplot)
##################################################

#####################
# EEG (1/f offset)
#####################

# Create line plot (EEG data)
eegplot_meanall = ggplot() +
  # Mean observed rt (actual data)
  geom_line(data = plotdat_ss, aes(y = plot_offset, x = trialnr), color = "#0072B2", size = 1.5) +
  # Posterior predicted PC (mean prediction)
  geom_line(data = plotdat_pc, aes(y = trial_eeg, x = trialnr), color = "#D55E00", size = 1.5) +
  # Posterior predicted rt (mean prediction)
  geom_line(data = plotdat_vp, aes(y = trial_eeg, x = trialnr), color = "#009E73", size = 1.5) +
  # Posterior predicted VP (mean prediction)
  geom_line(data = plotdat_ss, aes(y = trial_eeg, x = trialnr), color = "#CC79A7", size = 1.5) +
  labs(title = "PPC for all models (EEG)",
       x = "Trial Number",
       y = "1/f offset") +
  theme_classic(base_size = 24)+
  theme(
    plot.title = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "gray")
  )

eegplot1 = paste0("~\\EEGvar\\Visualization\\eegplot_meanall_", regions, ".png")
ggsave(plot = eegplot_meanall, eegplot1, dpi=300, height = 8, width = 12, bg = "white")

# Create line plot (RT model generated)
eegplot_MGall = ggplot() +
  # Posterior predicted PC (mean prediction)
  geom_line(data = plotdat_pc, aes(y = trial_eeg, x = trialnr), color = "#D55E00", size = 1.5) +
  # Posterior predicted rt (mean prediction)
  geom_line(data = plotdat_vp, aes(y = trial_eeg, x = trialnr), color = "#009E73", size = 1.5) +
  # Posterior predicted VP (mean prediction)
  geom_line(data = plotdat_ss, aes(y = trial_eeg, x = trialnr), color = "#CC79A7", size = 1.5) +
  labs(title = "PPC for all models (EEG)",
       x = "Trial Number",
       y = "1/f offset") +
  theme_classic(base_size = 24)+
  theme(
    plot.title = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "gray")
  )

eegplot2 = paste0("~\\EEGvar\\Visualization\\eegplot_MGall_", regions, ".png")
ggsave(plot = eegplot_MGall, eegplot2, dpi=300, height = 8, width = 12, bg = "white")

#####################
# Response times
#####################

# Create line plot (RT data)
rtplot_meanall = ggplot() +
  # Mean observed rt (actual data)
  geom_line(data = plotdat_ssrt, aes(y = plot_rt, x = trialnr), color = "#0072B2", size = 1.5) +
  # Posterior predicted PC (mean prediction)
  geom_line(data = plotdat_pcrt, aes(y = trial_rt, x = trialnr), color = "#D55E00", size = 1.5) +
  # Posterior predicted rt (mean prediction)
  geom_line(data = plotdat_vprt, aes(y = trial_rt, x = trialnr), color = "#009E73", size = 1.5) +
  # Posterior predicted VP (mean prediction)
  geom_line(data = plotdat_ssrt, aes(y = trial_rt, x = trialnr), color = "#CC79A7", size = 1.5) +
  labs(title = "PPC for all models (RT)",
       x = "Trial Number",
       y = "Response Time") +
  theme_classic(base_size = 24)+
  theme(
    plot.title = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "gray")
  )

# Create line plot (RT model generated)
rtplot_MGall = ggplot() +
  # Posterior predicted PC (mean prediction)
  geom_line(data = plotdat_pcrt, aes(y = trial_rt, x = trialnr), color = "#D55E00", size = 1.5) +
  # Posterior predicted rt (mean prediction)
  geom_line(data = plotdat_vprt, aes(y = trial_rt, x = trialnr), color = "#009E73", size = 1.5) +
  # Posterior predicted VP (mean prediction)
  geom_line(data = plotdat_ssrt, aes(y = trial_rt, x = trialnr), color = "#CC79A7", size = 1.5) +
  labs(title = "PPC for all models (RT)",
       x = "Trial Number",
       y = "Response Time") +
  theme_classic(base_size = 24)+
  theme(
    plot.title = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "gray")
  )

# plot names
rtplot1 = paste0("~\\EEGvar\\Visualization\\rtplot_meanall_", regions, ".png")
rtplot2 = paste0("~\\EEGvar\\Visualization\\rtplot_MGall_", regions, ".png")
# save plots
ggsave(plot = rtplot_meanall, rtplot1, dpi=300, height = 8, width = 12, bg = "black")
ggsave(plot = rtplot_MGall, rtplot2, dpi=300, height = 8, width = 12, bg = "black")

###########################################################
# Save model fits
###########################################################

# Model names
fit1name = paste0("~\\EEGvar\\Modeling\\pop_coding", regions, ".rds")
fit2name = paste0("~\\EEGvar\\Modeling\\var_precision_", regions, ".rds")
fit3name = paste0("~\\EEGvar\\Modeling\\seq_sensitive_", regions, ".rds")
# Save PC model
# saveRDS(fit1, file = fit1name)
# saveRDS(fit2, file = fit2name)
# saveRDS(fit3, file = fit3name)
# rm(fit1,fit2,fit3)
