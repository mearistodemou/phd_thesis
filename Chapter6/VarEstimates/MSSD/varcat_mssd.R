##########################################################
# Aggregate MSSD (whole brain & regional)
# Author: Michael E. Aristodemou
####################      Appendix      ##################
#1. Load dataframe & packages
#2. Create mean estimate of MSSD
#3. Create dataframe with MSSD across all conditions
#4. Region-wise data frame (cortical)
#5. Region-wise data frame (subcortical)
##########################################################

#####################
#1. Load dataframe
#####################

library(psych)
library(dplyr)
library(tidyr)
library(lavaan)

# Load voxel data frame
df1 = read.csv("~\\CatVar\\Data\\MSSD\\mssd_ses-drug1.csv")
df2 = read.csv("~\\CatVar\\Data\\MSSD\\mssd_ses-drug2.csv")
df3 = read.csv("~\\CatVar\\Data\\MSSD\\mssd_ses-drug3.csv")
# load region data frame (cortical)
df1_re = read.csv("~\\CatVar\\Data\\MSSD\\mssd_ses-drug1_regions.csv")
df2_re = read.csv("~\\CatVar\\Data\\MSSD\\mssd_ses-drug2_regions.csv")
df3_re = read.csv("~\\CatVar\\Data\\MSSD\\mssd_ses-drug3_regions.csv")
# load region data frame (subcortical)
df1_sub = read.csv("~\\CatVar\\Data\\MSSD\\mssd_ses-drug1_subreg.csv")
df2_sub = read.csv("~\\CatVar\\Data\\MSSD\\mssd_ses-drug2_subreg.csv")
df3_sub = read.csv("~\\CatVar\\Data\\MSSD\\mssd_ses-drug3_subreg.csv")

##########################################
#2. Create mean estimate of MSSD
##########################################

#########
# Drug1
#########

{
# Estimate mean MSSD per individual
df1 = df1 %>% group_by(Participant) %>%
      mutate(mssdBOLD = mean(MSSD, na.rm = TRUE)) %>%
      mutate(vox_count = row_number())

dfMSSD1 = df1[df1$vox_count == '1',] # create dataframe with mean MSSD
hist(dfMSSD1$mssdBOLD) # plot MSSD per person

# Extract the last number from the Participant string
dfMSSD1$subject <- as.numeric(sub("sub-0*", "", dfMSSD1$Participant))
dfMSSD1$session = 1

#########
# Drug 2
#########

# Estimate mean MSSD per individual
df2 = df2 %>% group_by(Participant) %>%
  mutate(mssdBOLD = mean(MSSD, na.rm = TRUE)) %>%
  mutate(vox_count = row_number())

dfMSSD2 = df2[df2$vox_count == '1',] # create dataframe with mean MSSD
hist(dfMSSD2$mssdBOLD) # plot MSSD per person

# Extract the last number from the Participant string
dfMSSD2$subject <- as.numeric(sub("sub-0*", "", dfMSSD2$Participant))
dfMSSD2$session = 2

#########
# Drug 3
#########

# Estimate mean MSSD per individual
df3 = df3 %>% group_by(Participant) %>%
  mutate(mssdBOLD = mean(MSSD, na.rm = TRUE)) %>%
  mutate(vox_count = row_number())

dfMSSD3 = df3[df3$vox_count == '1',] # create dataframe with mean MSSD
hist(dfMSSD3$mssdBOLD) # plot MSSD per person

# Extract the last number from the Participant string
dfMSSD3$subject <- as.numeric(sub("sub-0*", "", dfMSSD3$Participant))
dfMSSD3$session = 3

########################################################
#3. Create dataframe with MSSD across all conditions
########################################################

merged_df <- merge(dfMSSD1, dfMSSD2, by = c("subject"), all = TRUE)
merged_df <- merge(merged_df, dfMSSD3, by = c("subject"), all = TRUE)

# only keep variables of interest
fmriMID <- merged_df %>% select(subject, 
                                session, session.x, session.y, 
                                mssdBOLD, mssdBOLD.x, mssdBOLD.y)
fmriMID = 
  fmriMID %>%
  mutate(session.x = ifelse(is.na(session.x), 1, session.x)) %>%
  mutate(session.y = ifelse(is.na(session.y), 2, session.y)) %>%
  mutate(session = ifelse(is.na(session), 3, session))

# create long dataset
fmriMID_long <- fmriMID %>%
  pivot_longer(cols = c(session, session.x, session.y), names_to = "session_type", values_to = "session") %>%
  pivot_longer(cols = c(mssdBOLD, mssdBOLD.x, mssdBOLD.y), names_to = "mssdBOLD_type", values_to = "mssdBOLD") %>%
  filter(
    (session_type == "session" & mssdBOLD_type == "mssdBOLD") |
      (session_type == "session.x" & mssdBOLD_type == "mssdBOLD.x") |
      (session_type == "session.y" & mssdBOLD_type == "mssdBOLD.y")
  ) %>%
  select(subject, session, mssdBOLD)

# Load MID_df data
load("~\\CatVar\\Data\\MID_df.RData")

MID_df1 = MID_df[MID_df$trial == "1",]
MID_df1 = MID_df1 %>% select(subject, session, drug)


fmriMID = merge(MID_df1, fmriMID_long, by = c("subject", "session"), all = TRUE)

# save to data frame
write.csv(fmriMID, file = "~\\CatVar\\Data\\MSSD\\fmriMSSD.csv", row.names = FALSE)

}

##########################################################
#4. Region-wise data frame (cortical)
##########################################################

# Session number
df1_re$session = 1
df2_re$session = 2
df3_re$session = 3
# Participant id as number
df1_re$Participant <- as.numeric(sub("sub-0*", "", df1_re$Participant))
df2_re$Participant <- as.numeric(sub("sub-0*", "", df2_re$Participant))
df3_re$Participant <- as.numeric(sub("sub-0*", "", df3_re$Participant))

# merge dataframes
merged_df = merge(df1_re, df2_re, by = c("Participant", "Region"), all = TRUE)
merged_df = merge(merged_df, df3_re, by = c("Participant", "Region"), all = TRUE)

# only keep variables of interest
regionMID <- merged_df %>% select(Participant, session, session.x,
                                  session.y, MSSD, MSSD.x, MSSD.y, Region)
regionMID = 
  regionMID %>%
  mutate(session.x = ifelse(is.na(session.x), 1, session.x)) %>%
  mutate(session.y = ifelse(is.na(session.y), 2, session.y)) %>%
  mutate(session = ifelse(is.na(session), 3, session)) %>%
  mutate(subject = Participant)

# create long dataset
regionMID_long <- regionMID %>%
  pivot_longer(cols = c(session, session.x, session.y), names_to = "session_type", values_to = "session") %>%
  pivot_longer(cols = c(MSSD, MSSD.x, MSSD.y), names_to = "MSSD_type", values_to = "MSSD") %>%
  filter(
    (session_type == "session" & MSSD_type == "MSSD") |
      (session_type == "session.x" & MSSD_type == "MSSD.x") |
      (session_type == "session.y" & MSSD_type == "MSSD.y")
  ) %>%
  select(subject, session, MSSD, Region)

# Load MID_df data
load("~\\CatVar\\Data\\MID_df.RData")

MID_df1 = MID_df[MID_df$trial == "1",]
MID_df1 = MID_df1 %>% select(subject, session, drug)


regionMID = merge(MID_df1, regionMID_long, by = c("subject", "session"), all = TRUE)

# save to data frame
write.csv(regionMID, file = "~\\CatVar\\Data\\MSSD\\regionMSSD.csv", row.names = FALSE)


##########################################################
#5. Region-wise data frame (subcortical)
##########################################################

# Session number
df1_sub$session = 1
df2_sub$session = 2
df3_sub$session = 3
# Participant id as number
df1_sub$Participant <- as.numeric(sub("sub-0*", "", df1_sub$Participant))
df2_sub$Participant <- as.numeric(sub("sub-0*", "", df2_sub$Participant))
df3_sub$Participant <- as.numeric(sub("sub-0*", "", df3_sub$Participant))

# merge dataframes
merged_df = merge(df1_sub, df2_sub, by = c("Participant", "Region"), all = TRUE)
merged_df = merge(merged_df, df3_sub, by = c("Participant", "Region"), all = TRUE)

# only keep variables of interest
subMID <- merged_df %>% select(Participant, session, session.x,
                                  session.y, MSSD, MSSD.x, MSSD.y, Region)
subMID = 
  subMID %>%
  mutate(session.x = ifelse(is.na(session.x), 1, session.x)) %>%
  mutate(session.y = ifelse(is.na(session.y), 2, session.y)) %>%
  mutate(session = ifelse(is.na(session), 3, session)) %>%
  mutate(subject = Participant)

# create long dataset
subMID_long <- subMID %>%
  pivot_longer(cols = c(session, session.x, session.y), names_to = "session_type", values_to = "session") %>%
  pivot_longer(cols = c(MSSD, MSSD.x, MSSD.y), names_to = "MSSD_type", values_to = "MSSD") %>%
  filter(
    (session_type == "session" & MSSD_type == "MSSD") |
      (session_type == "session.x" & MSSD_type == "MSSD.x") |
      (session_type == "session.y" & MSSD_type == "MSSD.y")
  ) %>%
  select(subject, session, MSSD, Region)

# Load MID_df data
load("~\\CatVar\\Data\\MID_df.RData")

MID_df1 = MID_df[MID_df$trial == "1",]
MID_df1 = MID_df1 %>% select(subject, session, drug)


subMID = merge(MID_df1, subMID_long, by = c("subject", "session"), all = TRUE)

# save to data frame
write.csv(subMID, file = "~\\CatVar\\Data\\MSSD\\subMSSD.csv", row.names = FALSE)
