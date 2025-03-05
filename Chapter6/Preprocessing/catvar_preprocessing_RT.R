# Preprocessing for CatVar project (Chapter 6)
# Author: Michael E. Aristodemou
########################################################
# 1. Load datasets
# 2. Transform data into logRT
# 3. Split data into drug groups
# 4. Create MPLUS data
# 5. Preprocessing PET data
########################################################

########################################################
#1. Load datasets
########################################################
library(readr)
library(dplyr)

#Load dataframes
load("~\\CatVar\\Data\\MID_df.RData")
dat_pet <- read.csv("~\\CatVar\\Data\\roi_results_native.csv")
dat_demo <- read.table(file = "~\\CatVar\\Data\\all_proxies.tsv", header = TRUE, sep = "\t") %>%
  # one participant had an implausible estimate for sEBR (extreme outlier); treat this as NA
  dplyr::mutate(sEBR = ifelse(sEBR > 60, NA, sEBR)) %>%
  # code gender as a factor, and assign labels
  dplyr::mutate_at("gender", factor, labels = c("male", "female")) %>%
  dplyr::mutate_at("subject", factor)

######################################################
#2. Transform reaction times
######################################################

MID_df$logRT <- log(MID_df$RT_target)
hist(MID_df$RT_target)
hist(MID_df$logRT)

###################################
#3. Split data into drug groups
###################################
MID_df <- MID_df[, -which(names( MID_df) == "reward")]
MID_df$response <- ifelse(MID_df$response == "hit", 1, 0)

#3.1 Data subsetting according to drug used in trials
df_PBO <- MID_df[which(MID_df$drug == "PBO"),] #placebo trials
df_MPH <- MID_df[which(MID_df$drug == "MPH"),] #methylphenidate trials
df_SUL <- MID_df[which(MID_df$drug == "SUL"),] #sulpridite trials

# save files for R
write.csv(df_MPH, "~\\CatVar\\Data\\df_MPH.csv", row.names = FALSE)
write.csv(df_SUL, "~\\CatVar\\Data\\df_SUL.csv", row.names = FALSE)
write.csv(df_PBO, "~\\CatVar\\Data\\df_PBO.csv", row.names = FALSE)

####################################
#4. Create MPLUS data
####################################
#4.1. Transform NAs to 999 for MPLUS
df_MPH[is.na(df_MPH)] <- 999
df_SUL[is.na(df_SUL)] <- 999
df_PBO[is.na(df_PBO)] <- 999

#4.2 Save as .csv
write.csv(df_MPH, "~\\CatVar\\Mplus\\df_MPH_MPLUS.csv", row.names = FALSE)
write.csv(df_SUL, "~\\CatVar\\Mplus\\df_SUL_MPLUS.csv", row.names = FALSE)
write.csv(df_PBO, "~\\CatVar\\Mplus\\df_PBO_MPLUS.csv", row.names = FALSE)

############################################
#5. Preprocessing PET data
############################################
# load the PET Ki values
PET_ROI <- dat_pet %>%
  dplyr::select(subject, ROI, Ki = meanRoiContrastEstimate) %>%
  dplyr::filter(ROI %in% c("wholeStriatum", "wholePutamen", "wholeCaudate", "VS")) %>%
  dplyr::mutate(ROI = dplyr::recode_factor(
    ROI, 
    wholeStriatum = "striatum",
    wholePutamen = "putamen",
    wholeCaudate = "caudate",
    VS = "nacc")) %>%
  # make data frame wide with respect to the ROI
  tidyr::pivot_wider(id_cols = "subject", names_from = "ROI", values_from = "Ki", names_prefix = "Ki_")

write.csv(PET_ROI, "~\\CatVar\\PET_ROI.csv", row.names = FALSE)

