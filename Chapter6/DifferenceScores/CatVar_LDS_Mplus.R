######################################
# CatVar LDS models (Chapter 6)
# Using Mplus DSEM estimates
# Author: Michael E. Aristodemou
######################################

#####################################################
# 1. Install and load packages
# 2. Import factor scores from mplus
# 3. Create dataset with covariates
# 4. Preprocessing of data for LDS models 
# 5. Preprocessing of neural data
# 6. Specify bivariate LDS model MPH vs PBO
# 7. Specify bivariate LDS model SUL vs PBO
# 8. Visualize change scores under MPH and SUL
# 9. Rainclouds for neural variability change
# 10. Change score RTV x Ki_striatum covariance
# 11. Change score NEURAL x Ki_striatum covariance
# 12. Regional analyses  (cortical)
# 13. FDR correction for multiple comparisons
# 14. Extract standardized change for brain plots
# 15. Save regional change scores to csv
# 16. Bivariate LCS drift and RTV
#####################################################

########################################
#1. Install and load packages
########################################
{
  #1.1 Install packages
  list.of.packages <- c("MplusAutomation", "dplyr", "lavaan", "robustHD", "psych")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  #1.2 Load packages
  {
    library(MplusAutomation)
    library(dplyr)
    library(lavaan)
    library(robustHD)
    library(psych)
  }
}

load(file = "~\\CatVar\\catvarWD.RData")

##############################################
#2. Import factor scores from mplus file
##############################################

  #2.1 MPH factor scores
  #Rename factor scores
  {
    b2_factors <- read.table("~\\CatVar\\Mplus\\DSEM\\MPH\\20K\\DSEM_factors_MPH20K.dat", comment.char="")
    
    b2_factors <-
      b2_factors %>% 
      rename(
        
        LOGRT                            = V1,
        TRIAL                            = V2,
        LOGRT1                          = V3,
        PHI_Mean                         = V4,
        PHI_Median                       = V5,
        PHI_SD           = V6,
        PHI_LOWER                   = V7,
        PHI_HIGHER                  = V8,
        LOGV_Mean                        = V9,
        LOGV_Median                      = V10,
        LOGV_SD          = V11,
        LOGV_LOWER                  = V12,
        LOGV_HIGHER                 = V13,
        TREND_Mean                       = V14,
        TREND_Median                     = V15,
        TREND_SD         = V16,
        TREND_LOWER                 = V17,
        TREND_HIGHER                = V18,
        B_LOGRT_Mean                     = V19,
        B_LOGRT_Median                   = V20,
        B_LOGRT_SD       = V21,
        B_LOGRT_LOWER               = V22,
        B_LOGRT_HIGHER              = V23,
        SUBJECT                          = V24,
        TRUETIME                        = V25,
        TIMEPOINT                       = V26
        
      )
    
    #Make wide version of block 2 factors
    b2_wide <-b2_factors[which(b2_factors$TIMEPOINT == "1"),]
    
  }
  
  #2.2. Import factor scores from mplus file

    #2.2.1 MPH factor scores
    #Rename factor scores
    {
      b3_factors <- read.table("~\\CatVar\\Mplus\\DSEM\\PBO\\20K\\DSEM_factors_PBO20K.dat", comment.char="")
      
      b3_factors <-
        b3_factors %>% 
        rename(
          
          LOGRT                            = V1,
          TRIAL                            = V2,
          LOGRT1                          = V3,
          PHI_Mean                         = V4,
          PHI_Median                       = V5,
          PHI_SD           = V6,
          PHI_LOWER                   = V7,
          PHI_HIGHER                  = V8,
          LOGV_Mean                        = V9,
          LOGV_Median                      = V10,
          LOGV_SD          = V11,
          LOGV_LOWER                  = V12,
          LOGV_HIGHER                 = V13,
          TREND_Mean                       = V14,
          TREND_Median                     = V15,
          TREND_SD         = V16,
          TREND_LOWER                 = V17,
          TREND_HIGHER                = V18,
          B_LOGRT_Mean                     = V19,
          B_LOGRT_Median                   = V20,
          B_LOGRT_SD       = V21,
          B_LOGRT_LOWER               = V22,
          B_LOGRT_HIGHER              = V23,
          SUBJECT                          = V24,
          TRUETIME                        = V25,
          TIMEPOINT                       = V26
          
        )
      
      #Make wide version of block 2 factors
      b3_wide <-b3_factors[which(b3_factors$TIMEPOINT == "1"),]
      
    }
  
    #2.3 Import factor scores from mplus file
      #2.3.1 SUL factor scores
      #Rename factor scores
      {
        b4_factors <- read.table("~\\CatVar\\Mplus\\DSEM\\SUL\\20K\\DSEM_factors_SUL20K.dat", comment.char="")
        
        b4_factors <-
          b4_factors %>% 
          rename(
            
            LOGRT                            = V1,
            TRIAL                            = V2,
            LOGRT1                          = V3,
            PHI_Mean                         = V4,
            PHI_Median                       = V5,
            PHI_SD           = V6,
            PHI_LOWER                   = V7,
            PHI_HIGHER                  = V8,
            LOGV_Mean                        = V9,
            LOGV_Median                      = V10,
            LOGV_SD          = V11,
            LOGV_LOWER                  = V12,
            LOGV_HIGHER                 = V13,
            TREND_Mean                       = V14,
            TREND_Median                     = V15,
            TREND_SD         = V16,
            TREND_LOWER                 = V17,
            TREND_HIGHER                = V18,
            B_LOGRT_Mean                     = V19,
            B_LOGRT_Median                   = V20,
            B_LOGRT_SD       = V21,
            B_LOGRT_LOWER               = V22,
            B_LOGRT_HIGHER              = V23,
            SUBJECT                          = V24,
            TRUETIME                        = V25,
            TIMEPOINT                       = V26
            
          )
        
        #Make wide version of block 4 factors
        b4_wide <-b4_factors[which(b4_factors$TIMEPOINT == "1"),]
        
      }

######################################################################################
# 3. Create dataset with covariates
######################################################################################

# 3.1 Create dataset with PET covariates
  PET_ROI <- read.csv("~\\CatVar\\PET_ROI.csv")
  PET_ROI$SUBJECT <- PET_ROI$subject
  
# 3.2 Create dataset with ADHD covariates
  dat_demo <- read.table(file = "~\\CatVar\\Data\\all_proxies.tsv", header = TRUE, sep = "\t") %>%
    # one participant had an implausible estimate for sEBR (extreme outlier); treat this as NA
    dplyr::mutate(sEBR = ifelse(sEBR > 60, NA, sEBR)) %>%
    # code gender as a factor, and assign labels
    dplyr::mutate_at("gender", factor, labels = c("male", "female")) %>%
    dplyr::mutate_at("subject", factor)
  ADHD_col <- c("subject", "ADHD_adult", "ADHD_child")
  ADHD_df <- dat_demo[ADHD_col]
  ADHD_df$SUBJECT <- ADHD_df$subject
  
##########################################################################
# 4. Preprocessing of data for LDS models
##########################################################################

  # 4.1 Block comparison: PBO vs MPH

  {
    #b2 preprocessing for merge
    b2_wide$LOGV_Median_MPH <- b2_wide$LOGV_Median
    b2_wide$MEANRT_MPH<-b2_wide$B_LOGRT_Median
    my_vars <- c("SUBJECT", "LOGV_Median_MPH", "MEANRT_MPH")
    b2_merge <- b2_wide[my_vars]
    
    #b3 preprocessing for merge
    b3_wide$LOGV_Median_PBO <- b3_wide$LOGV_Median
    b3_wide$MEANRT_PBO<-b3_wide$B_LOGRT_Median
    my_vars <- c("SUBJECT", "LOGV_Median_PBO", "MEANRT_PBO")
    b3_merge <- b3_wide[my_vars]
    
    #Merge b2 & b3 by CLIENTID
    PBO_vs_MPH <- merge(b2_merge,b3_merge,by="SUBJECT")
    PBO_vs_MPH <- merge(PBO_vs_MPH, PET_ROI, by="SUBJECT")
    PBO_vs_MPH <- merge(PBO_vs_MPH, ADHD_df, by="SUBJECT")
    
    #Exponentiate
    #b2_vs_b3$b2_std <- exp(b2_vs_b3$LOGV_Median_b2)
    #b2_vs_b3$b3_std <- exp(b2_vs_b3$LOGV_Median_b3)
    
    #Rename
    PBO_vs_MPH$b2_PBO <- PBO_vs_MPH$LOGV_Median_PBO
    PBO_vs_MPH$b3_MPH <- PBO_vs_MPH$LOGV_Median_MPH
    PBO_vs_MPH$mean_PBO <- PBO_vs_MPH$MEANRT_PBO
    PBO_vs_MPH$mean_MPH <- PBO_vs_MPH$MEANRT_MPH
    
    
    #Rescale PET estimates
    PBO_vs_MPH$NACC<- scale(PBO_vs_MPH$Ki_nacc)
    PBO_vs_MPH$caudate<- scale(PBO_vs_MPH$Ki_caudate)
    PBO_vs_MPH$putamen<- scale(PBO_vs_MPH$Ki_putamen)
    PBO_vs_MPH$striatum<- scale(PBO_vs_MPH$Ki_striatum)
    
    #Rescale ADHD estimates
    PBO_vs_MPH$ADHD_ad <- scale(PBO_vs_MPH$ADHD_adult)
    PBO_vs_MPH$ADHD_cl <- scale(PBO_vs_MPH$ADHD_child)
  }
  
  
# 4.2 Block comparison: PBO vs SUL
  {

    #b3 preprocessing for merge
    b4_wide$LOGV_Median_SUL <- b4_wide$LOGV_Median
    b4_wide$MEANRT_SUL <- b4_wide$B_LOGRT_Median
    my_vars <- c("SUBJECT", "LOGV_Median_SUL", "MEANRT_SUL")
    b4_merge <- b4_wide[my_vars]
    
    #Merge b2 & b4 by CLIENTID
    PBO_vs_SUL <- merge(b3_merge,b4_merge,by="SUBJECT")
    PBO_vs_SUL <- merge(PBO_vs_SUL, PET_ROI, by="SUBJECT")
    PBO_vs_SUL <- merge(PBO_vs_SUL, ADHD_df, by="SUBJECT")
    
    #Exponentiate
    #b2_vs_b4$b2_std <- exp(b2_vs_b4$LOGV_Mean_b2)
    #b2_vs_b4$b4_std <- exp(b2_vs_b4$LOGV_Mean_b4)
    
    #Rename
    PBO_vs_SUL$b2_PBO <- PBO_vs_SUL$LOGV_Median_PBO
    PBO_vs_SUL$b4_SUL <- PBO_vs_SUL$LOGV_Median_SUL
    PBO_vs_SUL$mean_PBO <- PBO_vs_SUL$MEANRT_PBO
    PBO_vs_SUL$mean_SUL <- PBO_vs_SUL$MEANRT_SUL
    
    #Rescale PET estimates
    PBO_vs_SUL$NACC<- scale(PBO_vs_SUL$Ki_nacc)
    PBO_vs_SUL$caudate<- scale(PBO_vs_SUL$Ki_caudate)
    PBO_vs_SUL$putamen<- scale(PBO_vs_SUL$Ki_putamen)
    PBO_vs_SUL$striatum<- scale(PBO_vs_SUL$Ki_striatum)
    
    #Rescale ADHD estimates
    PBO_vs_SUL$ADHD_ad <- scale(PBO_vs_SUL$ADHD_adult)
    PBO_vs_SUL$ADHD_cl <- scale(PBO_vs_SUL$ADHD_child)
  }

#####################################
# 5. Preprocessing of neural data
#####################################
  
fmriMID = read.csv("~\\CatVar\\Data\\MSSD\\fmriMSSD.csv")
fmriWIDE = fmriMID %>% select(subject, drug, mssdBOLD)
fmriWIDE = tidyr::pivot_wider(fmriWIDE, names_from = drug, values_from = mssdBOLD)

# 5.1. Combine into 1 dataset for brain-behavior analyses
PBO_vs_MPH$subject = PBO_vs_MPH$SUBJECT
mixDATA = merge(fmriWIDE, PBO_vs_MPH, by = "subject", all = TRUE)
SULmerge = PBO_vs_SUL %>% select(SUBJECT, mean_SUL, b4_SUL) %>%
  rename(subject = SUBJECT)
mixDATA = merge(mixDATA, SULmerge, by = "subject", all = TRUE)
mixDATA$Ki_striatum = scale(mixDATA$Ki_striatum)
mixDATA$Ki_striatum_SQ = mixDATA$Ki_striatum * mixDATA$Ki_striatum

# 5.2. rename all
mixDATA <-
  mixDATA %>% 
  rename(neuSUL = SUL,
         neuMPH = MPH,
         neuPBO = PBO,
         rtvPBO = b2_PBO,
         rtvMPH = b3_MPH,
         rtvSUL = b4_SUL)

#######################################################################
# 6. Specify bivariate LDS model MPH vs PBO
#######################################################################

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

#######################################################################
# 7. Specify bivariate LDS model SUL vs PBO
#######################################################################

LCS_SULmix <- 
  '
#########################################
# NEURAL MODEL
#########################################

#Specify change score as latent variable
changeNEU =~ 1*neuSUL
#Perfectly regress b3 on b2
neuSUL ~ 1*neuPBO

changeNEU~~neuPBO #do not specify self-feedback

changeNEU ~ 1  #Change factor intercept

#Item intercepts
neuPBO ~ 1
neuSUL ~ 0*1 #equality constraint

changeNEU ~~ changeNEU  #Change factor variance

neuPBO ~~ neuPBO  
neuSUL ~~ 0*neuSUL


########################################
# BEHAVIORAL MODEL
########################################

#Specify change score as latent variable
changeRTV =~ 1*rtvSUL
#Perfectly regress b3 on b2
rtvSUL ~ 1*rtvPBO

changeRTV~~rtvPBO #do not specify self-feedback

changeRTV ~ 1  #Change factor intercept

#Item intercepts
rtvPBO ~ 1
rtvSUL ~ 0*1 #equality constraint

changeRTV ~~ changeRTV  #Change factor variance

rtvPBO ~~ rtvPBO  
rtvSUL ~~ 0*rtvSUL

##########################################
# COVARIATE MODEL
##########################################

# covariance between change scores
changeRTV ~~ changeNEU

# Baseline covariance
rtvPBO ~~ neuPBO

#ADD covariates for change
changeRTV ~ Ki_striatum + Ki_striatum_SQ
changeNEU ~ Ki_striatum + Ki_striatum_SQ

#COVARIATE intercepts
Ki_striatum ~ 1
Ki_striatum_SQ ~ 1

#COVARIATE variances
Ki_striatum ~~ Ki_striatum
Ki_striatum_SQ ~~ Ki_striatum_SQ

# Covary baseline Ki with baseline RTV and boldMSSD
Ki_striatum ~~ rtvPBO
Ki_striatum ~~ neuPBO

'
#fit model
fitLCS_SUL <- growth(LCS_SULmix, data=mixDATA, estimator='mlr', fixed.x=FALSE, missing='fiml')
summary(fitLCS_SUL, fit.measures=TRUE, standardized=TRUE)

#############################################################
# 8. Visualize change scores under SUL
#############################################################

# 8.1. MPH histograms
# get model predicted values
changemat_MPH<-lavPredict(fitLCS_MPH, type = "lv")

# histogram of MPH changes for neuVAR
changeMPH <- as.data.frame(changemat_MPH)
p1<-ggplot(changeMPH, aes(x=changeRTV)) + 
  geom_histogram(color="black", fill="#d19c2f", bins= 30)+
  theme_classic(base_size = 30)+
  labs(x = "RTV difference under MPH", y = "Frequency")
p1

ggsave(p1, file ='~\\CatVar\\Visualization\\MPH_dRTV_ch6.png', dpi = 300, height = 8, width = 10, bg= 'white')

# histogram of MPH changes for neuVAR
changeMPH <- as.data.frame(changemat_MPH)
p2<-ggplot(changeMPH, aes(x=changeNEU)) + 
  geom_histogram(color="black", fill="#d19c2f", bins= 30)+
  theme_classic(base_size = 30)+
  labs(x = "BOLD MSSD difference under MPH", y = "Frequency")
p2

ggsave(p2, file ='~\\CatVar\\Visualization\\MPH_dNEUVAR_ch6.png', dpi = 300, height = 8, width = 10, bg= 'white')


# 8.2. SUL histograms
# get model predicted values
changemat_SUL<-lavPredict(fitLCS_SUL, type = "lv")

# Histogram of SUL-driven change in boldMSSD
changeSUL <- as.data.frame(changemat_SUL)
p1<-ggplot(changeSUL, aes(x=changeRTV)) + 
  geom_histogram(color="black", fill="#49997c", bins= 30)+
  theme_classic(base_size = 30)+
  labs(x = "RTV difference under SUL", y = "Frequency")
p1

ggsave(p1, file ='~\\CatVar\\Visualization\\SUL_dRTV_ch6.png', dpi = 300, height = 8, width = 10, bg= 'white')

# Histogram of SUL-driven change in boldMSSD
changeSUL <- as.data.frame(changemat_SUL)
p2<-ggplot(changeSUL, aes(x=changeNEU)) + 
  geom_histogram(color="black", fill="#49997c", bins= 30)+
  theme_classic(base_size = 30)+
  labs(x = "BOLD MSSD difference under SUL", y = "Frequency")
p2

ggsave(p2, file ='~\\CatVar\\Visualization\\SUL_dNEUVAR_ch6.png', dpi = 300, height = 8, width = 10, bg= 'white')


##################################################################
# 9. Rainclouds for neural variability change
##################################################################

if (!require(remotes)) {
  install.packages("remotes")
}
remotes::install_github('njudd/ggrain')
library(ggrain)

#######################################
# 9.1. Rainclouds RTV
#######################################

######################################
#MPH versus placebo
######################################

#Pivot dataframe into long format
my_vars <- c("SUBJECT", "b2_PBO", "b3_MPH")
MPH_long <- PBO_vs_MPH[my_vars] #subset
MPH_long <- MPH_long %>% 
  rename(
    b_PBO = b2_PBO,
    b_MPH = b3_MPH
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
p2<-ggplot(MPH_long, aes(x = Drug, y = RTV, fill = Drug)) +
  geom_rain(alpha = .5, rain.side = 'f1x1', id.long.var = "SUBJECT",
            violin.args = list(color = NA, alpha = .7)) +
  theme_classic(base_size = 30) +
  #scale_x_discrete(limits = rev(levels(MPH_long$Drug)))+
  scale_fill_manual(values=c("#1ebecd", "#d19c2f")) +
  scale_color_manual(values=c("#1ebecd", "#d19c2f")) +
  ggsignif::geom_signif(
    comparisons = list(c("MPH", "PBO")),
    map_signif_level = TRUE) +
  guides(fill = 'none', color = 'none')

ggsave(p2, file ='~\\CatVar\\Visualization\\MPH_RTV_ch6.png', dpi = 300, height = 8, width = 10, bg= 'white')

######################################
#SUL versus placebo
######################################

#Pivot dataframe into long format
my_vars <- c("SUBJECT", "b2_PBO", "b4_SUL")
SUL_long <- PBO_vs_SUL[my_vars] #subset
SUL_long <- SUL_long %>% 
  rename(
    b_PBO = b2_PBO,
    b_SUL = b4_SUL
  )

#Pivot longer
SUL_long<- SUL_long %>%
  tidyr::pivot_longer(
    cols = starts_with("b"),
    names_to = "Drug",
    names_prefix = "b_",
    values_to = "RTV",
    values_drop_na = TRUE
  )

#Reverse factor order
SUL_long$Drug <- forcats::fct_rev(as.factor(SUL_long$Drug))

#Flanked rainclouds with jittered slopes
p3<-ggplot(SUL_long, aes(x = Drug, y = RTV, fill = Drug)) +
  geom_rain(alpha = .5, rain.side = 'f1x1', id.long.var = "SUBJECT",
            violin.args = list(color = NA, alpha = .7)) +
  theme_classic(base_size = 30) +
  #scale_x_discrete(limits = rev(levels(MPH_long$Drug)))+
  scale_fill_manual(values=c("#1ebecd", "#49997c")) +
  scale_color_manual(values=c("#1ebecd", "#49997c")) +
  ggsignif::geom_signif(
    comparisons = list(c("SUL", "PBO")),
    map_signif_level = TRUE) +
  guides(fill = 'none', color = 'none')
p3

ggsave(p3, file ='~\\CatVar\\Visualization\\SUL_RTV_ch6.png', dpi = 300, height = 8, width = 10, bg= 'white')

######################################
# 9.2. Rainclouds boldMSSD
######################################

#################################
# Methylphenidate
#################################

# Pivot dataframe into long format
my_vars <- c("subject", "neuPBO", "neuMPH")
MPH_long <- mixDATA[my_vars] #subset
MPH_long <- MPH_long %>% 
  rename(
    b_PBO = neuPBO,
    b_MPH = neuMPH
  )

# Pivot longer
MPH_long<- MPH_long %>%
  tidyr::pivot_longer(
    cols = starts_with("b"),
    names_to = "Drug",
    names_prefix = "b_",
    values_to = "RTV",
    values_drop_na = TRUE
  )

# Reverse factor order
MPH_long$Drug <- forcats::fct_rev(as.factor(MPH_long$Drug))

# Flanked rainclouds with jittered slopes
p2<-ggplot(MPH_long, aes(x = Drug, y = RTV, fill = Drug)) +
  geom_rain(alpha = .5, rain.side = 'f1x1', id.long.var = "subject",
            violin.args = list(color = NA, alpha = .7)) +
  theme_classic(base_size = 30) +
  ylab("BOLD MSSD") + xlab("Drug") +
  #scale_x_discrete(limits = rev(levels(MPH_long$Drug)))+
  scale_fill_manual(values=c("#1ebecd", "#d19c2f")) +
  scale_color_manual(values=c("#1ebecd", "#d19c2f")) +
  ggsignif::geom_signif(
    comparisons = list(c("MPH", "PBO")),
    map_signif_level = TRUE) +
  guides(fill = 'none', color = 'none')

ggsave(p2, file ='~\\CatVar\\Visualization\\MPH_NEUVAR_ch6.png', dpi = 300, height = 8, width = 10, bg= 'white')

######################################
# 9.2. SUL versus placebo NEURAL
######################################

# Pivot dataframe into long format
my_vars <- c("subject", "neuPBO", "neuSUL")
SUL_long <- mixDATA[my_vars] #subset
SUL_long <- SUL_long %>% 
  rename(
    b_PBO = neuPBO,
    b_SUL = neuSUL
  )

# Pivot longer
SUL_long<- SUL_long %>%
  tidyr::pivot_longer(
    cols = starts_with("b"),
    names_to = "Drug",
    names_prefix = "b_",
    values_to = "RTV",
    values_drop_na = TRUE
  )

# Reverse factor order
SUL_long$Drug <- forcats::fct_rev(as.factor(SUL_long$Drug))

# Flanked rainclouds with jittered slopes
p2<-ggplot(SUL_long, aes(x = Drug, y = RTV, fill = Drug)) +
  geom_rain(alpha = .5, rain.side = 'f1x1', id.long.var = "subject",
            violin.args = list(color = NA, alpha = .7)) +
  theme_classic(base_size = 30) +
  ylab("BOLD MSSD") + xlab("Drug") +
  #scale_x_discrete(limits = rev(levels(MPH_long$Drug)))+
  scale_fill_manual(values=c("#49997c", "#1ebecd")) +
  scale_color_manual(values=c("#49997c", "#1ebecd")) +
  ggsignif::geom_signif(
    comparisons = list(c("PBO", "SUL")),
    map_signif_level = TRUE) +
  guides(fill = 'none', color = 'none')

ggsave(p2, file ='~\\CatVar\\Visualization\\SUL_neuVAR_ch6.png', dpi = 300, height = 8, width = 10, bg= 'white')


########################################################################
# 10. Change score RTV x Ki_striatum covariance
########################################################################

#################################
# 10.1. Methylphenidate
#################################

#Create hacky diff scores
PBO_vs_MPH$diffMPH <- PBO_vs_MPH$b3_MPH - PBO_vs_MPH$b2_PBO
#Calculate sanity check correlation
cor(PBO_vs_MPH$diffMPH, PBO_vs_MPH$Ki_striatum)

#Plot
p5<-ggplot(PBO_vs_MPH, aes(x = Ki_striatum, y = diffMPH)) +
  geom_point(color="#d19c2f")+
  stat_smooth(method = "lm", color="#d19c2f")+
  theme_classic(base_size = 30)+
  labs(x = "Striatal [18F]-DOPA PET", y = "Change in RTV under MPH")
p5

ggsave(p5, file ='~\\CatVar\\Visualization\\MPH_PET_dRTV_ch6.png', dpi = 300, height = 8, width = 10, bg= 'white')

##################################
# 10.2. Sulpiride
##################################

#Create hacky diff scores
PBO_vs_SUL$diffSUL <- PBO_vs_SUL$b4_SUL - PBO_vs_SUL$b2_PBO
#Calculate sanity check correlation
cor(PBO_vs_SUL$diffSUL, PBO_vs_SUL$Ki_striatum)

#Plot
p6<-ggplot(PBO_vs_SUL, aes(x = Ki_striatum, y = diffSUL)) +
  geom_point(color="#49997c")+
  stat_smooth(method = "lm", color="#49997c")+
  theme_classic(base_size = 30)+
  labs(x = "Striatal [18F]-DOPA PET", y = "Change in RTV under SUL")
p6

ggsave(p6, file ='~\\CatVar\\Visualization\\SUL_PET_dRTV_ch6.png', dpi = 300, height = 8, width = 10, bg= 'white')

########################################################################
# 11. Change score NEURAL x Ki_striatum covariance
########################################################################

#################################
# 11.1 Methylphenidate
#################################

#Create hacky diff scores
mixDATA$diffMPH <- mixDATA$neuMPH - mixDATA$neuPBO
#Calculate sanity check correlation
cor(mixDATA$diffMPH, mixDATA$Ki_striatum)

#Plot
p5<-ggplot(mixDATA, aes(x = Ki_striatum, y = diffMPH)) +
  geom_point(color="#d19c2f")+
  stat_smooth(method = "lm", color="#d19c2f")+
  theme_classic(base_size = 30)+
  labs(x = "Striatal [18F]-DOPA PET", y = "Change in BOLD MSSD under MPH")
p5

ggsave(p5, file ='~\\CatVar\\Visualization\\MPH_PET_dNEU_ch6.png', dpi = 300, height = 8, width = 10, bg= 'white')

##################################
# 11.2 Sulpridite
##################################

#Create hacky diff scores
mixDATA$diffSUL <- mixDATA$neuSUL - mixDATA$neuPBO
#Calculate sanity check correlation
cor(mixDATA$diffSUL, mixDATA$Ki_striatum)

#Plot
p6<-ggplot(mixDATA, aes(x = Ki_striatum, y = diffSUL)) +
  geom_point(color="#49997c")+
  stat_smooth(method = "lm", color="#49997c")+
  theme_classic(base_size = 30)+
  labs(x = "Striatal [18F]-DOPA PET", y = "Change in BOLD MSSD under SUL")
p6

ggsave(p6, file ='~\\CatVar\\Visualization\\SUL_PET_dNEUVAR_ch6.png', dpi = 300, height = 8, width = 10, bg= 'white')

################################################################
# 12. Regional analyses  (cortical)
################################################################

#####################################
# 12.1. Preprocessing
#####################################

regionMID = read.csv("~\\CatVar\\Data\\MSSD\\regionMSSD.csv")
regionWIDE = regionMID %>% select(subject, drug, MSSD, Region)
regionWIDE = tidyr::pivot_wider(regionWIDE, names_from = drug, values_from = MSSD)
regionWIDE$Region = as.factor(regionWIDE$Region)

PBO_vs_MPH$subject = PBO_vs_MPH$SUBJECT
regionMIX = merge(regionWIDE, PBO_vs_MPH, by = "subject", all = TRUE)
SULmerge = PBO_vs_SUL %>% select(SUBJECT, mean_SUL, b4_SUL) %>%
  rename(subject = SUBJECT)
regionMIX = merge(regionMIX, SULmerge, by = "subject", all = TRUE)
regionMIX$Ki_striatum = scale(regionMIX$Ki_striatum)
regionMIX$Ki_striatum_SQ = regionMIX$Ki_striatum * regionMIX$Ki_striatum

# rename all
regionMIX <-
  regionMIX %>% 
  rename(neuSUL = SUL,
         neuMPH = MPH,
         neuPBO = PBO,
         rtvPBO = b2_PBO,
         rtvMPH = b3_MPH,
         rtvSUL = b4_SUL)

regionMIX = regionMIX[!is.na(regionMIX$Region),]


####################################################################
# 12.2. Run LSD per region (cortical)
####################################################################

#####################################
# Cortical regions
#####################################

# Methylphenidate
fitLCS_MPH = list()
pvalchange_MPH = list()
pvalcov_MPH = list()
estMPH_change = list()
p1 = list()

# Sulpiride
fitLCS_SUL = list()
pvalchange_SUL = list()
pvalcov_SUL = list()
estSUL_change = list()
p2 = list()

# Define regions to loop across
region_id = levels(regionMIX$Region)

for(i in region_id){
  
#fit model METHYLPHENIDATE
regionMIXsub = subset(regionMIX, Region == i)
fitLCS_MPH[[i]] <- growth(LCS_MPHmix, data=regionMIXsub, estimator='mlr', fixed.x=FALSE, missing='fiml')
p1[[i]] = parameterestimates(fitLCS_MPH[[i]], standardized = TRUE)
pvalchange_MPH[[i]] = p1[[i]]$pvalue[4]
pvalcov_MPH[[i]] = p1[[i]]$pvalue[19]
estMPH_change[[i]] = p1[[i]]$std.all[4]

#fit model SULPIRIDE
regionMIXsub = subset(regionMIX, Region == i)
fitLCS_SUL[[i]] <- growth(LCS_SULmix, data=regionMIXsub, estimator='mlr', fixed.x=FALSE, missing='fiml')
p2[[i]] = parameterestimates(fitLCS_SUL[[i]], standardized = TRUE)
pvalchange_SUL[[i]] = p2[[i]]$pvalue[4]
pvalcov_SUL[[i]] = p2[[i]]$pvalue[19]
estSUL_change[[i]] = p2[[i]]$std.all[4]

}

################################################################
# 12.3 Regional analyses  (subcortical)
################################################################

#####################################
# 12.4 Preprocessing
#####################################

subMID = read.csv("~\\CatVar\\Data\\MSSD\\subMSSD.csv")
subWIDE = subMID %>% select(subject, drug, MSSD, Region)
subWIDE = tidyr::pivot_wider(subWIDE, names_from = drug, values_from = MSSD)
subWIDE$Region = as.factor(subWIDE$Region)

PBO_vs_MPH$subject = PBO_vs_MPH$SUBJECT
subMIX = merge(subWIDE, PBO_vs_MPH, by = "subject", all = TRUE)
SULmerge = PBO_vs_SUL %>% select(SUBJECT, mean_SUL, b4_SUL) %>%
  rename(subject = SUBJECT)
subMIX = merge(subMIX, SULmerge, by = "subject", all = TRUE)
subMIX$Ki_striatum = scale(subMIX$Ki_striatum)
subMIX$Ki_striatum_SQ = subMIX$Ki_striatum * subMIX$Ki_striatum

# rename all
subMIX <-
  subMIX %>% 
  rename(neuSUL = SUL,
         neuMPH = MPH,
         neuPBO = PBO,
         rtvPBO = b2_PBO,
         rtvMPH = b3_MPH,
         rtvSUL = b4_SUL)

subMIX = subMIX[!is.na(subMIX$Region),]


####################################################################
# 12.5. Run LSD per region (subcortical)
####################################################################

#####################################
# Subcortical regions
#####################################

# Methylphenidate
fitLCS_MPHsub = list()
pvalchange_MPHsub = list()
pvalcov_MPHsub = list()
estMPH_changesub = list()

# Sulpiride
fitLCS_SULsub = list()
pvalchange_SULsub = list()
pvalcov_SULsub = list()
estSUL_changesub = list()

# Define regions to loop across
region_id = levels(subMIX$Region)

for(i in region_id){
  
  #fit model METHYLPHENIDATE
  subMIXsub = subset(subMIX, Region == i)
  fitLCS_MPHsub[[i]] <- growth(LCS_MPHmix, data=subMIXsub, estimator='mlr', fixed.x=FALSE, missing='fiml')
  p1[[i]] = parameterestimates(fitLCS_MPHsub[[i]], standardized = TRUE)
  pvalchange_MPHsub[[i]] = p1[[i]]$pvalue[4]
  pvalcov_MPHsub[[i]] = p1[[i]]$pvalue[19]
  estMPH_changesub[[i]] = p1[[i]]$std.all[4]
  
  #fit model SULPIRIDE
  subMIXsub = subset(subMIX, Region == i)
  fitLCS_SULsub[[i]] <- growth(LCS_SULmix, data=subMIXsub, estimator='mlr', fixed.x=FALSE, missing='fiml')
  p2[[i]] = parameterestimates(fitLCS_SULsub[[i]], standardized = TRUE)
  pvalchange_SULsub[[i]] = p2[[i]]$pvalue[4]
  pvalcov_SULsub[[i]] = p2[[i]]$pvalue[19]
  estSUL_changesub[[i]] = p2[[i]]$std.all[4]
  
}

###############################################
# 13. FDR correction for multiple comparisons
###############################################

# methylphenidate FDR
pvalMPH = unname(unlist(pvalchange_MPH)) # cortical
pvalMPHsub = unname(unlist(pvalchange_MPHsub)) # subcortical
pvalMPHcomb = c(pvalMPH, pvalMPHsub)
FDR_mph = p.adjust(pvalMPHcomb, method = "fdr", n = length(pvalMPHcomb))

# sulpiride FDR
pvalSUL = unname(unlist(pvalchange_SUL)) # cortical
pvalSULsub = unname(unlist(pvalchange_SULsub)) # subcortical
pvalSULcomb = c(pvalSUL, pvalSULsub)
FDR_sul = p.adjust(pvalSULcomb, method = "fdr", n = length(pvalSULcomb))

##################################################
# 14. Extract standardized change for brain plots
##################################################

# Cortical
############
# MPH change effect size
estMPH_change_df <- data.frame(
  Region = names(estMPH_change),
  Value = unlist(estMPH_change)
)
# MPH change pvalue
pvalchange_MPH_df = data.frame(
  Region = names(pvalchange_MPH),
  Value = head(FDR_mph, 48)
)
#SUL change efect size
estSUL_change_df <- data.frame(
  Region = names(estSUL_change),
  Value = unlist(estSUL_change)
)
# SUL change pvalue
pvalchange_SUL_df = data.frame(
  Region = names(pvalchange_SUL),
  Value = head(FDR_sul, 48)
)

# Subcortical
##############
# MPH change effect size
estMPH_changesub_df <- data.frame(
  Region = names(estMPH_changesub),
  Value = unlist(estMPH_changesub)
)
# MPH change pvalue
pvalchange_MPHsub_df = data.frame(
  Region = names(pvalchange_MPHsub),
  Value = tail(FDR_mph, 21)
)
# SUL change effect size
estSUL_changesub_df <- data.frame(
  Region = names(estSUL_changesub),
  Value = unlist(estSUL_changesub)
)
# SUL change pvalue
pvalchange_SULsub_df = data.frame(
  Region = names(pvalchange_SULsub),
  Value = tail(FDR_sul, 21)
)

##########################################
# 15. Save regional change scores to csv
##########################################

# std.all
# cortical MPH
write.csv(estMPH_change_df, "~\\CatVar\\Visualization\\MPH_delta_cort.csv", row.names = FALSE)  
# cortical SUL
write.csv(estSUL_change_df, "~\\CatVar\\Visualization\\SUL_delta_cort.csv", row.names = FALSE) 
# subcortical MPH
write.csv(estMPH_changesub_df, "~\\CatVar\\Visualization\\MPH_delta_subcort.csv", row.names = FALSE)  
# subcortical SUL
write.csv(estSUL_changesub_df, "~\\CatVar\\Visualization\\SUL_delta_subcort.csv", row.names = FALSE)

#pvalue
#cortical
write.csv(pvalchange_MPH_df, "~\\CatVar\\Visualization\\MPH_pval_cort.csv", row.names = FALSE)
write.csv(pvalchange_SUL_df, "~\\CatVar\\Visualization\\SUL_pval_cort.csv", row.names = FALSE)
#subcortical
write.csv(pvalchange_MPHsub_df, "~\\CatVar\\Visualization\\MPH_pval_subcort.csv", row.names = FALSE)
write.csv(pvalchange_SULsub_df, "~\\CatVar\\Visualization\\SUL_pval_subcort.csv", row.names = FALSE)

############################################
# 16. Bivariate LCS drift and RTV
############################################

library(dplyr)
subdrift <- participant_parameter_estimates %>% select(subject, drift_mph_mean, drift_pbo_mean)
driftDATA = merge(mixDATA, subdrift, by = "subject", all = TRUE)
driftDATA = driftDATA %>% rename(DRIFTMPH = drift_mph_mean,
                   DRIFTPBO = drift_pbo_mean)

LCS_MPHdrift <- 
  '
#########################################
# NEURAL MODEL
#########################################

#Specify change score as latent variable
changeDRIFT =~ 1*DRIFTMPH
#Perfectly regress b3 on b2
DRIFTMPH ~ 1*DRIFTPBO

changeDRIFT~~DRIFTPBO #do not specify self-feedback

changeDRIFT ~ 1  #Change factor intercept

#Item intercepts
DRIFTPBO ~ 1
DRIFTMPH ~ 0*1 #equality constraint

changeDRIFT ~~ changeDRIFT  #Change factor variance

DRIFTPBO ~~ DRIFTPBO  
DRIFTMPH ~~ 0*DRIFTMPH


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

changeRTV ~~ changeDRIFT # covariance between change scores
rtvPBO ~~ DRIFTPBO # baseline covariance

#ADD covariates for change
#changeRTV ~ Ki_striatum + Ki_striatum_SQ
#changeDRIFT ~ Ki_striatum + Ki_striatum_SQ

#COVARIATE intercepts
#Ki_striatum ~ 1
#Ki_striatum_SQ ~ 1

#COVARIATE variances
#Ki_striatum ~~ Ki_striatum
#Ki_striatum_SQ ~~ Ki_striatum_SQ

#Ki_striatum ~~ rtvPBO
#Ki_striatum ~~ DRIFTPBO


'
#fit model
fitLCS_MPH <- growth(LCS_MPHdrift, data=driftDATA, estimator='mlr', fixed.x=FALSE, missing='fiml')
summary(fitLCS_MPH, fit.measures=TRUE, standardized=TRUE)

modificationindices(fitLCS_MPH, sort. = TRUE)

#Extract individual difference-scores
changemat_MPH<-lavPredict(fitLCS_MPH , type = "lv")
changeDF = as.data.frame(changemat_MPH)
changeDF$subject <- 1:nrow(changeDF)

# Save as RDS file
saveRDS(changeDF, file = "changeDF.rds")


library(ggplot2)

ggplot(changemat_MPH, aes(x = changeDRIFT, y = changeRTV)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Lineplot of Value Over Time",
       x = "dDRIFT",
       y = "dRTV")

ggplot(changemat_MPH, aes(x = changeRTV)) +
  geom_histogram()

