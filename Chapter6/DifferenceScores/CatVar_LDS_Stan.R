######################################
# CatVar LDS models
# Using Stan DSEM estimates
# Author: Michael E. Aristodemou
######################################

library(dplyr)
library(lavaan)
library(tidyr)

#####################################################
# 1. Install and load packages
# 2. Import posterior estimates from Stan
# 3. Bivariate LDS neuVAR and RTV
#####################################################

datMPH = readRDS("~\\CatVar\\MID_collab\\Data\\Draws\\draws_MPH_standard.rds")
datPBO = readRDS("~\\CatVar\\MID_collab\\Data\\Draws\\draws_PBO_standard.rds")

############################################
# 2. Summarize posterior draws
############################################

# 2.1 MPH data formatting
# Assuming your data frame is named 'dat'
datMPH_long <- datMPH %>%
  summarize(across(everything(), ~ mean(.x, na.rm = TRUE)))
# save gammas
datMPH_gamma = datMPH_long %>%
  select(starts_with("gamma"))

# Only keep deviations (u)
datMPH_us = datMPH_long %>%
  select(starts_with("u")) %>%
  # Pivot the one-row data frame into long format
  pivot_longer(
    cols = everything(),
    names_to = c("id", "parameter"),
    names_pattern = "u\\.(\\d+)\\.(\\d+)",
    values_to = "mean_value"
  )

# Create wide data for analyses
datMPH_wide <- datMPH_us %>%
  pivot_wider(
    names_from = parameter,
    values_from = mean_value
  ) %>%
  rename(
    meanMPH =  "1",
    trendMPH="2",
    phiMPH="3",
    psiMPH="4"
  ) %>%
  mutate(meanMPH = meanMPH + datMPH_gamma$gamma.1,
         trendMPH = trendMPH + datMPH_gamma$gamma.2,
         phiMPH = phiMPH + datMPH_gamma$gamma.3,
         psiMPH = psiMPH + datMPH_gamma$gamma.4)


# 2.2 PBO data formatting
# Assuming your data frame is named 'dat'
datPBO_long <- datPBO %>%
  summarize(across(everything(), ~ mean(.x, na.rm = TRUE)))
# save gammas
datPBO_gamma = datPBO_long %>%
  select(starts_with("gamma"))

# Only keep us
datPBO_us = datPBO_long %>%
  select(starts_with("u")) %>%
  # Pivot the one-row data frame into long format
  pivot_longer(
    cols = everything(),
    names_to = c("id", "parameter"),
    names_pattern = "u\\.(\\d+)\\.(\\d+)",
    values_to = "mean_value"
  )

# Create wide data for analyses
datPBO_wide <- datPBO_us %>%
  pivot_wider(
    names_from = parameter,
    values_from = mean_value
  ) %>%
  rename(
    meanPBO =  "1",
    trendPBO ="2",
    phiPBO ="3",
    psiPBO ="4"
  ) %>%
  mutate(meanPBO = meanPBO + datPBO_gamma$gamma.1,
         trendPBO = trendPBO + datPBO_gamma$gamma.2,
         phiPBO = phiPBO + datPBO_gamma$gamma.3,
         psiPBO = psiPBO + datPBO_gamma$gamma.4)

# 2.3 Merge data frames
mixDATA <- left_join(datPBO_wide, datMPH_wide, by = "id") %>%
  select(starts_with(c("psi", "id"))) %>%
  rename(rtvMPH = psiMPH,
         rtvPBO = psiPBO) %>%
  mutate(subject = as.numeric(id))

# 2.4. Neural data preprocessing
fmriMID = read.csv("~\\CatVar\\Data\\MSSD\\fmriMSSD.csv")
fmriWIDE = fmriMID %>% select(subject, drug, mssdBOLD)
fmriWIDE = tidyr::pivot_wider(fmriWIDE, names_from = drug, values_from = mssdBOLD)
fmriWIDE = fmriWIDE %>% rename(neuSUL = SUL, neuMPH = MPH, neuPBO = PBO)
# 2.5. Merge data frames
mixDATA = left_join(fmriWIDE, mixDATA, by = "subject")

# 2.6. Add PET data
PET_ROI <- read.csv("~\\CatVar\\PET_ROI.csv")
mixDATA = left_join(PET_ROI, mixDATA, by = "subject")
mixDATA$Ki_striatum = scale(mixDATA$Ki_striatum)
mixDATA$Ki_striatum_SQ = mixDATA$Ki_striatum * mixDATA$Ki_striatum

###########################################
# 3. Bivariate LDS neuVAR and RTV
###########################################

# 3.1. MPH vs PBO model
LCS_MPHmix <- 
'
#########################################
# NEURAL MODEL
#########################################

#Specify change score as latent variable
changeNEU =~ 1*neuMPH
#Perfectly regress b3 on b2
neuMPH ~ 1*neuPBO

changeNEU~~neuPBO #do not specify self-feedback

changeNEU ~ 1  #Change factor intercept

#Item intercepts
neuPBO ~ 1
neuMPH ~ 0*1 #equality constraint

changeNEU ~~ changeNEU  #Change factor variance

neuPBO ~~ neuPBO  
neuMPH ~~ 0*neuMPH


########################################
# BEHAVIORAL MODEL
########################################

#Specify change score as latent variable
changeRTV =~ 1*rtvMPH
#Perfectly regress b3 on b2
rtvMPH ~ 1*rtvPBO

changeRTV~~rtvPBO #do not specify self-feedback

changeRTV ~ 1  #Change factor intercept

#Item intercepts
rtvPBO ~ 1
rtvMPH ~ 0*1 #equality constraint

changeRTV ~~ changeRTV  #Change factor variance

rtvPBO ~~ rtvPBO  
rtvMPH ~~ 0*rtvMPH

##########################################
# COVARIATE MODEL
##########################################

# covariance between change scores
changeRTV ~~ changeNEU

#ADD covariates for change
changeRTV ~ Ki_striatum + Ki_striatum_SQ
changeNEU ~ Ki_striatum + Ki_striatum_SQ

#COVARIATE intercepts
Ki_striatum ~ 1
Ki_striatum_SQ ~ 1

#COVARIATE variances
Ki_striatum ~~ Ki_striatum
Ki_striatum_SQ ~~ Ki_striatum_SQ

rtvPBO ~~ neuPBO
Ki_striatum ~~ rtvPBO
Ki_striatum ~~ neuPBO

'
#fit model
fitLCS_MPH <- growth(LCS_MPHmix, data=mixDATA, estimator='mlr', fixed.x=FALSE, missing='fiml')
summary(fitLCS_MPH, fit.measures=TRUE, standardized=TRUE)

modificationindices(fitLCS_MPH, sort. = TRUE)

#Extract individual difference-scores
changemat_SUL<-lavPredict(fitLCS_MPH, type = "lv")

#################################
# Visualization: neuVAR x RTV
#################################
library(ggplot2)

#################
# 1. Histograms
#################

# Change in RTV after MPH
p1<-ggplot(changemat_MPH, aes(x=changeRTV)) + 
  geom_histogram(color="black", fill="#d19c2f", bins= 30)+
  theme_classic(base_size = 26)+
  labs(x = "RTV difference under MPH", y = "Frequency")
p1

# Change in neuVAR after MPH
p2<-ggplot(changemat_MPH, aes(x=changesup)) + 
  geom_histogram(color="#d19c2f", fill="black", bins= 30)+
  theme_classic(base_size = 26)+
  labs(x = "Surprise difference under MPH", y = "Frequency")
p2

ggsave(p1, file ='~\\CatVar\\Chapter6\\Visualization\\hist_changeRTVsurp_MPH.png', 
       dpi = 300, height = 8, width = 10, bg= 'white')
ggsave(p2, file ='~\\CatVar\\Chapter6\\Visualization\\hist_changeSUP_MPH.png', 
       dpi = 300, height = 8, width = 10, bg= 'white')

######################
# 2. Regression plot
######################

# Regression plot for change-change correlation
p3 = ggplot(changemat_MPH, aes(x = changesup, y = changeRTV)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color="#d19c2f") +
  theme_classic(base_size = 18)+
  labs(title = "Negative dSurprise ~~ dRTV",
       x = "dSurprise",
       y = "dRTV")
p3

ggsave(p3, file ='~\\CatVar\\Chapter6\\Visualization\\regplot_dRTVxdSUP_MPH.png', 
       dpi = 300, height = 8, width = 10, bg= 'white')


########################
# 3. Rainclouds
########################

if (!require(remotes)) {
  install.packages("remotes")
}
remotes::install_github('njudd/ggrain')
library(ggrain)

# 3.1 RTV: MPH vs PBO
#Pivot dataframe into long format
my_vars <- c("id", "rtvPBO", "rtvMPH")
MPH_long <- surpDATA[my_vars] #subset
MPH_long <- MPH_long %>% 
  rename(
    b_PBO = rtvPBO,
    b_MPH = rtvMPH
  )

#Pivot longer
MPH_long<- MPH_long %>%
  tidyr::pivot_longer(
    cols = starts_with("b"),
    names_to = "Drug",
    names_prefix = "b_",
    values_to = "RTV",
    values_drop_na = TRUE
  )

#Reverse factor order
MPH_long$Drug <- forcats::fct_rev(as.factor(MPH_long$Drug))

#Flanked rainclouds with jittered slopes
p4<-ggplot(MPH_long, aes(x = Drug, y = RTV, fill = Drug)) +
  geom_rain(alpha = .5, rain.side = 'f1x1', id.long.var = "id",
            violin.args = list(color = NA, alpha = .7)) +
  theme_classic(base_size = 30) +
  #scale_x_discrete(limits = rev(levels(MPH_long$Drug)))+
  scale_fill_manual(values=c("#1ebecd", "#d19c2f")) +
  scale_color_manual(values=c("#1ebecd", "#d19c2f")) +
  ggsignif::geom_signif(
    comparisons = list(c("MPH", "PBO")),
    map_signif_level = TRUE) +
  guides(fill = 'none', color = 'none')

# save raincloud dRTV MPH vs PBO 
ggsave(p4, file ='~\\CatVar\\Chapter6\\Visualization\\rain_RTVsurp_MPH.png', 
       dpi = 300, height = 8, width = 10, bg= 'white')

# 3.2 DRFIT: MPH vs PBO
#Pivot dataframe into long format
my_vars <- c("id", "supPBO", "supMPH")
MPH_long <- surpDATA[my_vars] #subset
MPH_long <- MPH_long %>% 
  rename(
    b_PBO = supPBO,
    b_MPH = supMPH
  )

#Pivot longer
MPH_long<- MPH_long %>%
  tidyr::pivot_longer(
    cols = starts_with("b"),
    names_to = "Drug",
    names_prefix = "b_",
    values_to = "Surprise",
    values_drop_na = TRUE
  )

#Reverse factor order
MPH_long$Drug <- forcats::fct_rev(as.factor(MPH_long$Drug))

#Flanked rainclouds with jittered slopes
p5<-ggplot(MPH_long, aes(x = Drug, y = Surprise, fill = Drug)) +
  geom_rain(alpha = .5, rain.side = 'f1x1', id.long.var = "id",
            violin.args = list(color = NA, alpha = .7)) +
  theme_classic(base_size = 30) +
  #scale_x_discrete(limits = rev(levels(MPH_long$Drug)))+
  scale_fill_manual(values=c("#1ebecd", "#d19c2f")) +
  scale_color_manual(values=c("#1ebecd", "#d19c2f")) +
  ggsignif::geom_signif(
    comparisons = list(c("MPH", "PBO")),
    map_signif_level = TRUE) +
  guides(fill = 'none', color = 'none')

# save raincloud dRTV MPH vs PBO 
ggsave(p5, file ='~\\CatVar\\Chapter6\\Visualization\\rain_SUP_MPH.png', 
       dpi = 300, height = 8, width = 10, bg= 'white')
